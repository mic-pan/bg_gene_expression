# Option to run simulations, or load previous results.
run_benchmark = false

using BondGraphs, ModelingToolkit, VegaLite, DataFrames, DifferentialEquations, Sundials, CSV
using BenchmarkTools, Statistics
include("gene_regulation.jl")
include("vl_functions.jl")

sym2var(x) = (@variables $x(t))[1]
vec2df(x,y) = DataFrame(x=x,y=y)

# Define parameters
μA_nom = 20
A_nom = 5800000
γ_max = 1260
kf_eff = 4*γ_max # Multiplied by 4 to account for 4 ATP molecules per amino acid. Effective rate constant because it is for the reaction C_i ⇌ C_i+1, where A has been absorbed into parameters.
α = exp(1.25)
kt = 100*kf_eff

## Define functions
"""Set parameters for full translation model."""
function set_params_translation!(sys,n;μ_A=μA_nom)
    sys.A.e = μ_A
    sys.M.K = 1.0
    sys.M.q = 10.0
    sys.R.K = 1.0
    sys.R.q = 5000.0
    sys.C0.K = 1.0
    sys.C0.q = 1e-6

    sys.C1.K = α/(sys.M.K*sys.M.q)
    sys.C1.q = 1e-6
    sys.r1.r = kf_eff/exp(μA_nom)
    intermediates = Dict(
        i => getproperty(sys.elong,Symbol("X$(i-1)")) for i in 2:n-1
    )
    intermediates[1] = sys.C1
    intermediates[n] = getproperty(sys,Symbol("C$n"))

    for i in 2:n
        K_prev = intermediates[i-1].K
        c = intermediates[i]
        set_species_param!(c,:K,α*K_prev)
        c.q = 1e-6
        getproperty(sys.elong,Symbol("r$(i-1)")).r = kf_eff/K_prev/exp(μA_nom)
    end

    set_reaction_param!(sys.rd, :r, kt/get_species_param(intermediates[n],:K))
    sys.rb.r = 1e-3
    sys.P.e = -15

    return sys
end

"""Set parameters for elongation module"""
function set_params_simple_elongation!(B,elong,KC_base;μ_A=μA_nom)
    m = elong.n
    B.K = KC_base/m
    B.q = 1e-6*m

    elong.α = α
    elong.r = kf_eff/KC_base*α/exp(μA_nom)
end

"""Set parameters for simple translation model."""
function set_params_translation_simple!(sys;μ_A=μA_nom)
    n = sys.elong.n+1

    sys.M.K = 1.0
    sys.M.q = 10.0
    sys.R.K = 1.0
    sys.R.q = 5000.0
    
    sys.C0.K = 1.0
    sys.C0.q = 1e-6
    sys.r1.r = kf_eff/exp(μA_nom)
    sys.A.e = μ_A

    K_C1 = α/(sys.M.K*sys.M.q)
    set_params_simple_elongation!(sys.B0,sys.elong,K_C1;μ_A=μA_nom)
    
    sys.Cn.μ0 = n*log(α) - log(sys.M.K*sys.M.q)
    sys.Cn.q = 1e-6

    sys.rd.μr = log(kt) - sys.Cn.μ0
    sys.rb.r = 1e-3
    sys.P.e = -15

    return sys
end

"""Set parameters for translation with multiple lumped complexes."""
function set_params_translation_blocks!(sys;μ_A=μA_nom)
    sys.M.K = 1.0
    sys.M.q = 10.0
    sys.R.K = 1.0
    sys.R.q = 5000.0
    
    sys.C0.K = 1.0
    sys.C0.q = 1e-6
    sys.r1.r = kf_eff/exp(μA_nom)
    sys.A.e = μ_A

    n_el = sum(type.(nodes(sys.elong.bondgraph)) .== "re_el")
    K_base = α/(sys.M.K*sys.M.q)
    for i in 1:n_el
        lumped_complex = (i==1) ? sys.B0 : getproperty(sys.elong,Symbol("B$(i-1)"))
        re_el = getproperty(sys.elong,Symbol("el$i"))
        set_params_simple_elongation!(lumped_complex,re_el,K_base;μ_A=μA_nom)
        K_base *= α^re_el.n
    end
    
    n = sum(getproperty(sys.elong,Symbol("el$i")).n for i in 1:n_el)+1
    sys.Cn.μ0 = n*log(α) - log(sys.M.K*sys.M.q)
    sys.Cn.q = 1e-6

    sys.rd.μr = log(kt) - sys.Cn.μ0
    sys.rb.r = 1e-3
    sys.P.e = -15
    return sys
end

"""Plot comparison of dynamics between models."""
function plot_lines(df,xsymb,ysymb;xlabel="x",ylabel="y",width=120,height=90)
    return df |> @vlplot(
        :line, width=width, height=height,
        x={xsymb, axis={title=xlabel}}, 
        y={ysymb, axis={title=ylabel}},
        color={
            :model, 
            sort=["1 complex","2 complexes","5 complexes","Full model"],
            scale={range=["#f462a3", "#bb61a5", "#7d5da7", "#1657a8"]},
        },
        row=:n,
        column={:A_norm,title="A (10⁶ molecules)"},
        resolve={scale={x=:independent,y=:independent}},
        config={
            view={stroke=:black},
            font="Arial",
            legend={title=""},
            axis={
                title="",
                titleFontWeight=:normal,
                grid=false,
                tickSize=3,
                domainColor=:black,
                tickColor=:black,
            }
        }
    )
end

"""Swap Ce components to Se components for translation model."""
function swap_to_se!(sys)
    A = Component(:se,:A)
    swap!(sys,sys.A,A)
    P = Component(:se,:P)
    swap!(sys,sys.P,P)
    return sys
end

"""Simulate translation for different models."""
function simulate_translation(n,A)
    print("n=$n, A=$A\n")
    sys = Translation(n=n,open=false,n_lumps=0)
    swap_to_se!(sys)
    set_params_translation!(sys,n)
    
    sys_simple = Translation(n=n,open=false,n_lumps=1)
    swap_to_se!(sys_simple)
    set_params_translation_simple!(sys_simple)

    sys_blocks_2 = Translation(n=n,open=false,n_lumps=2)
    swap_to_se!(sys_blocks_2)
    set_params_translation_blocks!(sys_blocks_2)

    sys_blocks_5 = Translation(n=n,open=false,n_lumps=5)
    swap_to_se!(sys_blocks_5)
    set_params_translation_blocks!(sys_blocks_5)

    systems =  [sys,sys_simple,sys_blocks_2,sys_blocks_5]
    
    # Test effects at different energy levels
    μA = 20 .+ log(A/A_nom)
    for sys in systems
        sys.A.e = μA
    end

    tmax = round(5800*n/A,sigdigits=1)
    tspan = (0.0,tmax)
    options = Dict(:solver => CVODE_BDF(), :reltol => 1e-9, :abstol => 1e-12)

    sols = [simulate(sys,tspan;options...) for sys in systems]
    return sols
end

"""Convert ODE solution to DataFrame."""
function sol2df(x)
    n = x[1][1]
    A_norm = x[1][2]/1e6
    (sol,sol_simple,sol_blocks2,sol_blocks5) = x[2]

    @variables C0₊q(t) C1₊q(t) B0₊q(t) elong₊B1₊q(t) Cn₊q(t) R₊q(t) rd₊p1₊F(t) elong₊v(t) elong₊el2₊v(t)
    macroexpand(Main, :(@variables C0₊q(t)) )
    _syms = [Symbol("elong₊X$(i)₊q") for i in 1:n-2]
    syms = [_syms; :C0₊q; :C1₊q; Symbol("C$(n)₊q")]
    complex_vars = sym2var.(syms)
    total_complex = sum(complex_vars)

    _syms = [Symbol("elong₊B$(i)₊q") for i in 1:10]
    syms = [:B0₊q; _syms]
    lumped_complex_vars = sym2var.(syms)

    t_arr = range(0,sol.t[end],length=201)
    df_full = DataFrame(n=n,A_norm=A_norm,model=:"Full model",t=t_arr,
        v=sol(t_arr,idxs=rd₊p1₊F).u/1000,c=sol(t_arr,idxs=total_complex).u,r=sol(t_arr,idxs=R₊q).u)
    df_simple = DataFrame(n=n,A_norm=A_norm,model="1 complex",t=t_arr,
        v=sol_simple(t_arr,idxs=rd₊p1₊F/1000).u,c=sol_simple(t_arr,idxs=C0₊q+B0₊q+Cn₊q).u,r=sol_simple(t_arr,idxs=R₊q).u)
    df_blocks2 = DataFrame(n=n,A_norm=A_norm,model="2 complexes",t=t_arr,
        v=sol_blocks2(t_arr,idxs=rd₊p1₊F/1000).u,c=sol_blocks2(t_arr,idxs=C0₊q+B0₊q+elong₊B1₊q+Cn₊q).u,r=sol_blocks2(t_arr,idxs=R₊q).u)
    df_blocks5 = DataFrame(n=n,A_norm=A_norm,model="5 complexes",t=t_arr,
        v=sol_blocks5(t_arr,idxs=rd₊p1₊F/1000).u,c=sol_blocks5(t_arr,idxs=C0₊q+sum(lumped_complex_vars[1:5])+Cn₊q).u,r=sol_blocks5(t_arr,idxs=R₊q).u)

    return [df_full;df_simple;df_blocks2;df_blocks5]
end

## Run dynamic simulations
n_vec = [20,60,100]
A_vec = [100000,1000000,5800000]
sols = [((n,A),simulate_translation(n,A)) for n in n_vec for A in A_vec];

## Plot comparison figures
dfs = [sol2df(x) for x in sols]
df_all = vcat(dfs...)
fig_v = plot_lines(df_all,:t,:v;xlabel="Time (mins)",ylabel="Translation rate")
fig_c = plot_lines(df_all,:t,:c;xlabel="Time (mins)",ylabel="Total complex")
fig_r = plot_lines(df_all,:t,:r;xlabel="Time (mins)",ylabel="Free ribosome")

savevl(fig_v,"output/v_comparison")
savevl(fig_c,"output/c_comparison")
savevl(fig_r,"output/r_comparison")

## Benchmark run times for simple and complex models
n_vec = 20:20:400
full_model_times = BenchmarkTools.Trial[]
simple_model_times = BenchmarkTools.Trial[]
block2_times = BenchmarkTools.Trial[]
block5_times = BenchmarkTools.Trial[]

"""Run a benchmark of the ODE system."""
function benchmark_model(sys)
    tmax = 50.0
    tspan = (0.0,tmax)
    options = Dict(:reltol => 1e-9, :abstol => 1e-12)
    μA = 20
    sys.A.e = μA
    alg = CVODE_BDF()
    odesys = ODESystem(sys)
    prob = ODEProblem(odesys, [], tspan; reltol=1e-9, abstol=1e-12)
    return @benchmark solve($prob, $alg) seconds=60 samples=1000
end

"""Returns the mean and standard deviation of runtimes."""
function extract_times(times,scaling_factor=1e6)
    μ = [mean(x.times) for x in times]/scaling_factor
    σ = [std(x.times) for x in times]/scaling_factor
    return (μ,σ)
end

df_benchmark(n_vec,μ,σ,label) = DataFrame(
    n=n_vec,μ=μ,σ=σ,t1=μ-σ,t2=μ+σ,model=label
)

if run_benchmark
    # Run simulations
    for n in n_vec
        print("n=$n\n")
        sys = Translation(n=n,open=false,n_lumps=0)
        swap_to_se!(sys)
        set_params_translation!(sys,n)

        sys_simple = Translation(n=n,open=false,n_lumps=1)
        swap_to_se!(sys_simple)
        set_params_translation_simple!(sys_simple)

        sys_blocks_2 = Translation(n=n,open=false,n_lumps=2)
        swap_to_se!(sys_blocks_2)
        set_params_translation_blocks!(sys_blocks_2)
    
        sys_blocks_5 = Translation(n=n,open=false,n_lumps=5)
        swap_to_se!(sys_blocks_5)
        set_params_translation_blocks!(sys_blocks_5)
        
        models = [sys,sys_simple,sys_blocks_2,sys_blocks_5]
        (t1,t2,t3,t4) = benchmark_model.(models)

        push!(full_model_times,t1)
        push!(simple_model_times,t2)
        push!(block2_times,t3)
        push!(block5_times,t4)
    end

    # Extract run times
    scaling_factor = 1e6
    (μ_full_model,σ_full_model) = extract_times(full_model_times,scaling_factor)
    (μ_simple_model,σ_simple_model) = extract_times(simple_model_times,scaling_factor)
    (μ_blocks2,σ_blocks2) = extract_times(block2_times,scaling_factor)
    (μ_blocks5,σ_blocks5) = extract_times(block5_times,scaling_factor)

    # Compile results into DataFrame
    df_full = df_benchmark(n_vec,μ_full_model,σ_full_model,"Full model")
    df_simple = df_benchmark(n_vec,μ_simple_model,σ_simple_model,"1 complex")
    df_blocks2 = df_benchmark(n_vec,μ_blocks2,σ_blocks2,"2 complexes")
    df_blocks5 = df_benchmark(n_vec,μ_blocks5,σ_blocks5,"5 complexes")
    df = vcat(df_full,df_simple,df_blocks2,df_blocks5)

    # Save results
    CSV.write("output/benchmark.csv",df)
else
    # Load results
    df = DataFrame(CSV.File("output/benchmark.csv"))
end

# Plot comparison of runtimes between full and simple models
fig = df |> vl_config() + @vlplot(
    :line, width=180, height=135,
    x={:n, title="n"}, 
    y={
        :μ, 
        scale={type=:log}, 
        axis={
            title="Mean execution time (ms)",
            grid=true
        }
    },
    color={
        :model, 
        sort=["1 complex","2 complexes","5 complexes","Full model"],
        scale={range=["#f462a3", "#bb61a5", "#7d5da7", "#1657a8"]},
    }
)
savevl(fig,"output/benchmark")

## Simulate energy dissipation for a regular length protein
A_vec = [100000,1000000,5800000]
n = 1200
sys_simple = Translation(n=n,open=false,n_lumps=1)
swap_to_se!(sys_simple)
set_params_translation_simple!(sys_simple)
odesys = ODESystem(sys_simple)

@parameters A₊e P₊e
@variables rd₊p1₊F(t) rb₊p1₊E(t) rb₊p2₊E(t) elong₊p1₊E(t) elong₊p2₊E(t) elong₊p3₊E(t) rd₊p1₊E(t) rd₊p2₊E(t) r1₊p1₊E(t) r1₊p2₊E(t)
translation_efficiency(μA,μP,n) = μP/(n*μA)
μA_vec = @. 20 + log(A_vec/A_nom)
μ_folding = 20.0
μP = 1.25*n - μ_folding + log(1000) # Assuming protein amount of 1000 molecules
efficiency = translation_efficiency.(μA_vec,μP,n)

# Compile DataFrame with results
df = DataFrame(A=A_vec,efficiency=efficiency,v=NaN,μA=μA_vec,μP=μP,
    ΔG_binding=NaN,ΔG_elongation=NaN,ΔG_termination=NaN)
    tspan = (0.0,500.0)
for (i,μA) in enumerate(μA_vec)
    ps = [
        A₊e => μA,
        P₊e => μP
    ]
    prob = ODEProblem(odesys, [], tspan, ps)
    sol = solve(prob, AutoVern7(Rodas5()))
    df.v[i] = sol[rd₊p1₊F][end]
    df.ΔG_binding[i] = sol[rb₊p2₊E-rb₊p1₊E][end]
    df.ΔG_elongation[i] = sol[elong₊p3₊E-elong₊p1₊E-(n-1)*elong₊p2₊E + r1₊p2₊E-r1₊p1₊E][end]
    df.ΔG_termination[i] = sol[rd₊p2₊E-rd₊p1₊E][end]
end

# Modify DataFrame for presentation as a table
df_display = copy(df)
df_display.efficiency = round.(100*df_display.efficiency,sigdigits=3)
df_display.v = round.(df_display.v,digits=1)
df_display.μA = round.(df_display.μA,digits=1)
df_display.μP = round.(df_display.μP,digits=1)
df_display.ΔG_binding = round.(df_display.ΔG_binding,digits=1)
df_display.ΔG_elongation = round.(df_display.ΔG_elongation,digits=1)
df_display.ΔG_termination = round.(df_display.ΔG_termination,sigdigits=3)