using DataFrames, VegaLite, Polynomials, LinearAlgebra
include("vl_functions.jl")

# Define parameters
rb = 1e-2
KM = 1
KR = 1
dM = log(2)/2
dP = log(2)/4
r_Tc = 1.13e-2
R_bA = 366.4
h = 2
b = 100
KPh_RbI = b^h

"""Define the polynomial to be solved in finding the steady states"""
function ss_polynomial(R)
    a = rb*KM*KR*R*r_Tc*R_bA/(dM*dP)
    return Polynomial([-a*b^2,a^2+b^2,-2a,2,-a/b^2,1/b^2])
end

"""Solve the for the steady-state values of protein"""
function toggle_ss(R)
    p = roots(ss_polynomial(R))
    p_real = p[isreal.(p)]
    P_ss = [x.re for x in p_real]
    return P_ss
end

P1 = 10 .^ range(-1,stop=log10(4000),length=100)
P2_base = @. rb*KM*KR/(dM*dP) * r_Tc*R_bA/(1 + P1^h/KPh_RbI)

# Find the fixed points for different ribosome numbers
df = DataFrame(R=Float64[], P1=Float64[], P2=Float64[], nullcline=Symbol[])
df_fixed_points = DataFrame(R=Float64[], P1=Float64[], P2=Float64[], type=Symbol[])
R_vec = [1000,2000,5000]

for R in R_vec
    P2 = R*P2_base
    df_new1 = DataFrame(R=R, P1=P1, P2=P2, nullcline=:P2)
    append!(df, df_new1)
    df_new2 = DataFrame(R=R, P1=P2, P2=P1, nullcline=:P1)
    append!(df, df_new2)
    P_ss = toggle_ss(R)
    df_fixed_new = DataFrame(
        R=R, 
        P1=P_ss[[1,3,2]], 
        P2=P_ss[[3,1,2]], 
        type=[:stable, :stable, :unstable]
    )
    append!(df_fixed_points, df_fixed_new)
end
df_fixed_points.text = [
    "($(Int(round(P1))),$(Int(round(P2))))" 
    for (P1,P2) in zip(df_fixed_points.P1,df_fixed_points.P2)
]

# Plot the phase planes
function plot_phase_plane(df,df_fixed_points,R)
    df_subset = df[df.R .== R,:]
    df_fixed_subset= df_fixed_points[df_fixed_points.R .== R,:]
    df_stable_subset = df_fixed_subset[df_fixed_subset.type .== :stable,:]
    df_unstable_subset = df_fixed_subset[df_fixed_subset.type .== :unstable,:]

    return @vlplot() + @vlplot(
        :line, data = df_subset,
        title="R = $R",
        x="P1:q",
        y="P2:q",
        color={:nullcline, scale={range=["#1657a8", "#048059"]}},
    ) + @vlplot(
        {:circle, color="#6626bb", size=40, stroke=:white, strokeWidth=1, opacity=1}, 
        data = df_stable_subset,
        x="P1:q",
        y="P2:q"
    ) + @vlplot(
        {:circle, color="#e03e0e", size=40, stroke=:white, strokeWidth=1, opacity=1}, 
        data = df_unstable_subset,
        x="P1:q",
        y="P2:q",
        text=:type
    ) + @vlplot(
        {:text, color="#6626bb", dx=4, dy=-4, align=:left, baseline=:bottom},
        data = df_stable_subset,
        x="P1:q",
        y="P2:q",
        text=:text,
    ) + @vlplot(
        {:text, color="#e03e0e", dx=1, dy=-7, align=:left, baseline=:bottom},
        data = df_unstable_subset,
        x="P1:q",
        y="P2:q",
        text=:text,
    )
end

figs = [plot_phase_plane(df,df_fixed_points,R) for R in R_vec]
fig_phase_plane = @vlplot(config={font="Arial"}) + [figs[1] figs[2] figs[3]]
savevl(fig_phase_plane,"output/toggle_phase_plane")
