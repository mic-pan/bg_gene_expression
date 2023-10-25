load_results = false

using Peaks, DataFrames, Sundials, CSV, LinearAlgebra, Plots, VegaLite
include("gene_circuits.jl")
include("vl_functions.jl")
default(:fontfamily, "Helvetica")

cA_nom = 9.63e-3
A_nom = 5800000

## Run the gene expression model for different energy and ribosome levels
model = GeneExpression(open=false)
A = Component(:se,:A)
swap!(model,model.A,A)
set_params_gr!(model)
I = Component(:se,:I)
swap!(model,model.I,I)
μ0_P = model.P.μ0
I.e = μ0_P + log(1e-3) # low inhibitor concentration
sys = ODESystem(model)

@parameters A₊e
@variables Tl₊elong₊p1₊F(t)
tspan = (0.0,500.0)
options = Dict(
    :solver => CVODE_BDF(),
    :abstol => 1e-12,
    :reltol => 1e-9
)

A_nom = 5800000
A_vec = 10000:10000:A_nom
v_tl_vec = zeros(size(A_vec))

for (i,xA) in enumerate(A_vec)
    ps = [A₊e => 20 + log(xA/A_nom)]
    prob = ODEProblem(sys, [], tspan, ps; options...)
    sol = solve(prob,CVODE_BDF())
    v_tl_vec[i] = sol[Tl₊elong₊p1₊F][end]
end

v_tl_ref = v_tl_vec[end]
v_half_max = v_tl_ref/2
idx = findlast(v_tl_vec .< v_half_max)
t1 = A_vec[idx]
x1 = v_tl_vec[idx]
t2 = A_vec[idx+1]
x2 = v_tl_vec[idx+1]
K_A = (t2-t1)/(x2-x1)*(v_half_max-x1) + t1
K_A_norm = K_A/1000
K_A_conc = cA_nom *K_A/A_nom

df_A = DataFrame(A=A_vec/1000,v=v_tl_vec)
df_A_subset = subset(df_A, :A => ByRow(<=(200)))

fig_gr_A1 = df_A |> vl_config() + @vlplot(
    {:line, color="#1657a8"}, width=120, height=90,
    x={:A, title="A (10³ molecules)"}, 
    y={:v, title="Translation rate"}
)

fig_gr_A2 = df_A_subset |> vl_config() + @vlplot(
    {:line, color="#1657a8"}, width=120, height=90,
    x={:A, title="A (10³ molecules)"}, 
    y={:v, title="Translation rate"}
) + @vlplot(
    {:circle, color="#1657a8"},
    x={datum=K_A_norm},
    y={datum=v_half_max},
) + @vlplot(
    {:rule, color=:grey, strokeDash=[4,4]},
    x={datum=K_A_norm}
) + @vlplot(
    {:text, align=:left, dx=10},
    x={datum=K_A_norm},
    y={datum=v_half_max},
    text={datum="Kₘ = $(round(K_A_norm,digits=1))"}
)

@variables R₊q(t)
R_vec = 100:100:5000
v_tl_vec = zeros(size(R_vec))

for (i,R) in enumerate(R_vec)
    u0 = [R₊q => R]
    prob = ODEProblem(sys, u0, tspan, []; options...)
    sol = solve(prob,CVODE_BDF())
    v_tl_vec[i] = sol[Tl₊elong₊p1₊F][end]
end

df_R = DataFrame(R=R_vec,v=v_tl_vec)
fig_gr_R = df_R |> vl_config() + @vlplot(
    {:line, color="#1657a8"}, width=120, height=90,
    x={:R, title="Ribosomes"}, 
    y={:v, title="Translation rate"}
)

savevl(fig_gr_A1,"output/gr_A1")
savevl(fig_gr_A2,"output/gr_A2")
savevl(fig_gr_R,"output/gr_R")

## Toggle switch
toggle_results_path = "output/toggle_resources.csv"

model = ToggleSwitch()
set_params_toggle!(model)
A = Component(:se,:A)
swap!(model,model.A,A)
model.P1.q = 100.0
model.P2.q = 10.0
sys = ODESystem(model)

tspan = (0.0,500.0)

A_vec = 1000*(1000:50:6000)
R_vec = 1000:50:6000
@parameters A₊e
@variables G1₊mRNA₊q(t) G2₊mRNA₊q(t) G3₊mRNA₊q(t) P1₊q(t) P2₊q(t) P3₊q(t) R₊q(t)
if load_results == false
    df = DataFrame(A=Float64[], R=Float64[], S=Float64[], τ=Float64[])
    for xA in A_vec, R in R_vec
        u0 = [R₊q => R]
        ps = [A₊e => 20 + log(xA/A_nom)]
        prob = ODEProblem(sys, u0, tspan, ps; options...)
        sol = solve(prob,CVODE_BDF())

        P1 = sol[P1₊q][end]
        P2 = sol[P2₊q][end]
        sensitivity = P1/P2
        p_ss = [P1,P2]

        p_array = sol[[P1₊q,P2₊q]]
        dist_norm = [p./p_ss .- 1 for p in p_array]
        dist_vec = norm.(dist_norm)
        target_dist = 0.01
        idx = findlast(dist_vec .> target_dist)

        t1 = sol.t[idx]
        x1 = dist_vec[idx]
        t2 = sol.t[idx+1]
        x2 = dist_vec[idx+1]
        response_time = (t2-t1)/(x2-x1)*(target_dist-x1) + t1

        #print("A = $xA, R = $R , S = $sensitivity, τ = $response_time\n")
        df_new = DataFrame(A=xA/1e3,R=R,S=sensitivity,τ=response_time)
        append!(df,df_new)
    end
    CSV.write(toggle_results_path, df)
else
    df = DataFrame(CSV.File(toggle_results_path))
end

p1 = contour(
    R_vec,A_vec/1000,reshape(df.S,(length(R_vec),length(A_vec)))',
    color=:viridis, fill=true, lw=0, size=(300,300), aspect_ratio=1, 
    labelfontsize=10, titlefontsize=11, c=:haline
)
xlabel!("Ribosomes")
ylabel!("A (10³ molecules)")
xlims!(1000,6000)
ylims!(1000,6000)
title!("Bistability index")

p2 = contour(
    R_vec,A_vec/1000,reshape(df.τ,(length(R_vec),length(A_vec)))',
    color=:viridis, fill=true, lw=0, size=(300,300), aspect_ratio=1, 
    labelfontsize=10, titlefontsize=11, c=:haline
)
xlabel!("Ribosomes")
ylabel!("A (10³ molecules)")
xlims!(1000,6000)
ylims!(1000,6000)
title!("Switching time (mins)")

l = @layout [a b]
p = plot(p1,p2, layout=l, size=(600,300))
savefig(p,"output/toggle_resources.svg")
savefig(p,"output/toggle_resources.pdf")

## Repressilator
repressilator_results_path = "output/repressilator_resources.csv"

model = Repressilator()
set_params_repressilator!(model)
A = Component(:se,:A)
swap!(model,model.A,A)
sys = ODESystem(model)

tspan = (0.0,1500.0)

A_vec = 1000*(1000:50:6000)
R_vec = 1000:50:6000

if load_results == false
    df = DataFrame(A=Float64[], R=Float64[], amplitude=Float64[], period=Float64[])
    @variables G1₊mRNA₊q(t) G2₊mRNA₊q(t) G3₊mRNA₊q(t) P1₊q(t) P2₊q(t) P3₊q(t) R₊q(t)
    for xA in A_vec, R in R_vec
        A.e = 20 + log(xA/A_nom)
        model.R.q = R
        u0 = [R₊q => R]
        ps = [A₊e => 20 + log(xA/A_nom)]
        prob = ODEProblem(sys, u0, tspan, ps; options...)
        sol = solve(prob,CVODE_BDF())

        t_vec = 1000:0.1:1500
        P_vec = sol(t_vec,idxs=P1₊q)

        (idx_max,peaks) = findmaxima(P_vec)
        (idx_min,troughs) = findminima(P_vec)
        t_peaks = t_vec[idx_max]

        if length(t_peaks) > 2 && t_peaks[end-1] > 100.0
            period = t_peaks[end] - t_peaks[end-1]
            amplitude = peaks[end] - troughs[end]
            if amplitude < 1.0
                amplitude = NaN
                period = NaN
            end
        else
            period = NaN
            amplitude = NaN
        end

        print("A = $xA, R = $R, Amplitude = $amplitude, Period = $period \n")
        df_new = DataFrame(A=xA,R=R,amplitude=amplitude,period=period)
        append!(df,df_new)
    end
    CSV.write(repressilator_results_path, df)
else
    df = DataFrame(CSV.File(repressilator_results_path))
end

p1 = contour(
    R_vec,A_vec/1000,reshape(df.amplitude,(length(R_vec),length(A_vec)))',
    color=:viridis, fill=true, lw=0, size=(300,300), aspect_ratio=1, 
    labelfontsize=10, titlefontsize=11, c=:haline
)
xlabel!("Ribosomes")
ylabel!("A (10³ molecules)")
xlims!(1000,6000)
ylims!(1000,6000)
title!("Amplitude")

p2 = contour(
    R_vec,A_vec/1000,reshape(df.period,(length(R_vec),length(A_vec)))',
    color=:viridis, fill=true, lw=0, size=(300,300), aspect_ratio=1, 
    labelfontsize=10, titlefontsize=11, c=:haline
)
xlabel!("Ribosomes")
ylabel!("A (10³ molecules)")
xlims!(1000,6000)
ylims!(1000,6000)
title!("Period (mins)")

l = @layout [a b]
p = plot(p1,p2, layout=l, size=(600,300))
savefig(p,"output/repressilator_resources.svg")
savefig(p,"output/repressilator_resources.pdf")
