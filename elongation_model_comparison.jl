using BondGraphs, ModelingToolkit, VegaLite, DataFrames, DifferentialEquations
include("gene_regulation.jl")
include("vl_functions.jl")

## Define functions
"""Swaps dynamic components to source of effort components in an elongation model."""
function swap_to_se!(sys,m)
    A = Component(:se,:A)
    swap!(sys,sys.A,A)
    X0 = Component(:se,:X0)
    swap!(sys,sys.X0,X0)
    Xm = Component(:se,"X$m")
    swap!(sys,getproperty(sys,Symbol("X$m")),Xm)
    return sys
end

"""Set parameters for elongation model."""
function set_params_elongation!(sys,m)
    sys.X0.e = log(X0)

    for i in 1:m-1
        getproperty(sys,Symbol("X$i")).K = 1.0
        getproperty(sys,Symbol("X$i")).q = 0.1
    end
    for i in 1:m
        getproperty(sys,Symbol("r$i")).r = 1.0
    end

    getproperty(sys,Symbol("X$m")).e = log(Xm)
    
    return sys
end

"""Plot comparison of values between models"""
function plot_lines(df,xsym,ysym,csym;xlabel="",ylabel="")
    return df |> vl_config() + @vlplot(
        :line, width=180, height=135,
        x={xsym, axis={title = xlabel}}, 
        y={ysym, axis={title = ylabel}},
        color={csym, scale={range=["#1657a8", "#f462a3"]}},
        strokeDash={csym, scale={range=[[4,0],[4,4]]}}
    )
end

## Elongation parameters
m = 5
kf = 1.0
kr = 1.0
X0 = 1
Xm = 0.1

## Code dynamic model of elongation using bond graphs
sys = Elongation(m=m,open=false)
swap_to_se!(sys,m)
set_params_elongation!(sys,m)

A_vec = 0.01:0.01:2.0
tspan = (0.0,50.0)
sols = []

## Simulate dynamic model under different energy levels
for A in A_vec
    sys.A.e = log(A)
    sol = simulate(sys,tspan,solver=Rodas5())
    push!(sols,sol)
end

# Record elongation rate from simulations
@variables r5₊p1₊F(t) 
v_el_vec = [sol[r5₊p1₊F][end] for sol in sols]

## Plot comparison of elongation rates
A_norm = kf*A_vec/kr
v_analytic = @. kr*(A_norm^m*X0 - Xm)*(1-A_norm)/(1-A_norm^m)
i = findfirst(A_norm .== 1)
v_analytic[i] = kr*(X0 - Xm)/m
df1 = DataFrame(A=A_norm,v=v_el_vec,model=:Full)
df2 = DataFrame(A=A_norm,v=v_analytic,model=:Reduced)
df = [df1;df2]

fig_v_vs_A = plot_lines(df,:A,:v,:model,xlabel="Â",ylabel="Elongation rate")
savevl(fig_v_vs_A,"output/v_trans_comparison")

## Plot comparison of complex amounts
sym2var(x) = (@variables $x(t))[1]
syms = [Symbol("X$(i)₊q") for i in 1:m-1]
complex_vars = sym2var.(syms)

C = [[sol[c][end] for sol in sols] for c in complex_vars]

a = @. -(X0*A_norm^m - Xm)/(1 - A_norm^m)
b = @. (X0 - Xm)/(1 - A_norm^m)
C_analytic = [@. a + b*A_norm^i for i in 1:m-1]

dfs_full = [DataFrame(A = A_norm, c = c_vec, index = i, model=:Full) for (i,c_vec) in enumerate(C)]
dfs_reduced = [DataFrame(A = A_norm, c = c_vec, index = i, model=:Reduced) for (i,c_vec) in enumerate(C_analytic)]
df = vcat(dfs_full...,dfs_reduced...)

fig_c_vs_A = df |> vl_config() + @vlplot(
    :line, width=180, height=135,
    x={:A, axis={title = "Â"}}, 
    y={:c, axis={title = "Amount"}},
    detail = :index,
    color = {:model, scale={range=["#1657a8", "#f462a3"]}},
    strokeDash={:model, scale={range=[[4,0],[4,4]]}}
)
savevl(fig_c_vs_A,"output/c_comparison_elongation")

## Run simulations for different chain lengths
m_vec = 2:40
A = 1.2
v_el_vec = Float64[]
for m in m_vec
    sys = Elongation(m=m,open=false)
    swap_to_se!(sys,m)
    sys.A.e = log(A)
    set_params_elongation!(sys,m)

    tspan = (0.0,1000.0)
    sol = simulate(sys,tspan,solver=Rodas5())

    v = sym2var(Symbol("r$(m)₊p1₊F"))
    push!(v_el_vec, sol[v][end])
end

## Plot results for different chain lengths
A_norm = kf*A/kr
v_analytic = @. kr*(A_norm^m_vec*X0 - Xm)*(1-A_norm)/(1-A_norm^m_vec)
df1 = DataFrame(m=m_vec,v=v_el_vec,model=:Full)
df2 = DataFrame(m=m_vec,v=v_analytic,model=:Reduced)
df = [df1;df2]
fig_m_vs_A = plot_lines(df,:m,:v,:model,xlabel="m",ylabel="Elongation rate")
savevl(fig_m_vs_A,"output/m_comparison")

## Calculate analytical elongation rate for different chain lengths
m_vec = 2 .^ (1:6)
A_norm = 0:0.01:2
function trans_rate(A_norm,m,X0,Xm,kr)
    if A_norm == 1
        return kr*(X0-Xm)/m
    else
        return kr*(A_norm^m*X0 - Xm)*(1-A_norm)/(1-A_norm^m)
    end
end

function trans_rate_df(A_norm,m,X0,Xm,kr)
    v_analytic = trans_rate.(A_norm,m,X0,Xm,kr)
    df = DataFrame(A_norm=A_norm,v=v_analytic,m=m)
    return df
end

df_sims = vcat(collect(trans_rate_df(A_norm,m,X0,Xm,kr) for m in m_vec)...)

v_el_approx = [(A < 1) ? -kr*Xm*(1-A) : kr*X0*(A-1) for A in A_norm]
df_approx = DataFrame(A_norm=A_norm,v=v_el_approx,m=Inf)

## Plot elongation rate against energy for different chain lengths
fig_v_vs_n = df_sims |> vl_config() + @vlplot(
    :line, width=180, height=135,
    x={:A_norm, axis={title = "Â"}}, 
    y={:v,axis={title = "Elongation rate"}},
    color={"m:n",scale={range=["#aed0ff", "#91afe5", "#738fcc", "#5570b3", "#35539b", "#023883"]}}
)
savevl(fig_v_vs_n,"output/v_vs_n")