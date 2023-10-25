# Option to run simulations or load results from previous simulations
load_results = true

using Peaks, DataFrames, Sundials, LinearAlgebra, Distributions, Random, VegaLite, CSV
using JSON, Optim, Statistics
include("gene_circuits.jl")
include("vl_functions.jl")

# Basic plot options
fig_config = @vlplot(
    resolve={scale={color=:independent}},
    config={
        view={stroke=:transparent},
        font="Arial",
        legend={disable=true},
        axis={
            titleFontWeight=:normal,
            grid=false,
            tickSize=3,
            domainColor=:black,
            tickColor=:black,
        },
    }
)

"""Returns the amplitude and period from oscillatory differential equation solution."""
function analyse_oscillations(sol)
    (idx_max,peaks) = findmaxima(sol[P1₊q])
    (idx_min,troughs) = findminima(sol[P1₊q])
    t_peaks = sol.t[idx_max]

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
    return (amplitude,period)
end

"""Add dots corresponding to specific points to a distribution."""
quantile_plot(x) = @vlplot(
    {:point, opacity=1},
    data = {values=[
        {x=x[1]},
        {x=x[2]},
        {x=x[3]},
        {x=x[4]},
        {x=x[5]},
    ]},
    x="x:q",
    y={datum=0},
    color={"x:n",scale={scheme=:viridis}}
)

"""Plot timeseries of protein number for specific parameter values."""
function plot_timeseries(percentiles,sols)
    df_timeseries = DataFrame(t=Float64[], P=Float64[], Percentile=Float64[])
    t_plot = 0:1:250
    for (p,sol) in zip(percentiles,sols)
        df_new = DataFrame(t=t_plot, P=sol(t_plot,idxs=P1₊q).u, Percentile=p)
        append!(df_timeseries,df_new)
    end

    fig_timeseries = df_timeseries |> fig_config + @vlplot(
        :line, width=120,height=90,
        x = {:t, title="Time (mins)"},
        y = {:P, title="Protein 1"},
        color = {"Percentile:n", scale={scheme=:viridis}},
    )
    fig_timeseries.params["config"]["legend"]["disable"] = false
    return (fig_timeseries,df_timeseries)
end

## Define repressilator model
model = Repressilator()
set_params_repressilator!(model)
A = Component(:se,:A)
swap!(model,model.A,A)
A.e = 20
sys = ODESystem(model; simplify_eqs=true)

tspan = (0.0,500.0)
options = Dict(
    :solver => CVODE_BDF(),
    :abstol => 1e-9,
    :reltol => 1e-6
)

## Simulate effects of heterogeneity in protein degradation rates
Pdeg_results_path = "output/repressilator_distribution.csv"
Random.seed!(5688308)
@parameters G1₊P_deg₊deg₊μr G2₊P_deg₊deg₊μr G3₊P_deg₊deg₊μr
@variables G1₊mRNA₊q(t) G2₊mRNA₊q(t) G3₊mRNA₊q(t) P1₊q(t) P2₊q(t) P3₊q(t) R₊q(t)
function simulate_degradation(sys,μr)
    ps = [
        G1₊P_deg₊deg₊μr => μr,
        G2₊P_deg₊deg₊μr => μr,
        G3₊P_deg₊deg₊μr => μr
    ]
    prob = ODEProblem(sys, [], tspan, ps; options...)
    sol = solve(prob, CVODE_BDF())

    return sol
end

# Run dynamic simulations
if load_results == false
    μ_deg_mean = model.G1.P_deg.deg.μr
    μ_deg_std = log(2)/4 # To achieve 2-fold change within a 95% confidence interval
    d = Normal(μ_deg_mean,μ_deg_std)
    μr_deg_vec = rand(d, 10000)

    df = DataFrame(half_life=Float64[], amplitude=Float64[], period=Float64[])
     for μr in μr_deg_vec
        half_life = log(2)/exp(μr+model.P1.μ0)
        sol = simulate_degradation(sys,μr)
        (amplitude,period) = analyse_oscillations(sol)
        df_new = DataFrame(half_life=half_life,amplitude=amplitude,period=period)
        append!(df,df_new)
    end

    CSV.write(Pdeg_results_path, df)
else
    df = DataFrame(CSV.File(Pdeg_results_path))
end

# Extract parameters in 5, 25, 50, 75, 95th percentiles
percentiles = [5,25,50,75,95]
half_life_vals = [quantile(df.half_life,p/100) for p in percentiles]
μr_vec = @. log(log(2)/half_life_vals) - model.P1.μ0

# Run dynamic simulations for specific parameters
sols = [simulate_degradation(sys,μr) for μr in μr_vec];
array_oscillations = analyse_oscillations.(sols)
amplitudes = [a[1] for a in array_oscillations]
periods = [a[2] for a in array_oscillations]

# Plot distributions
fig_half_life = df |> @vlplot() +  @vlplot(
    {:area, color="#a373e3"}, width=120,height=90,
    transform=[
        {density="half_life",bandwidth=0.2}
    ],
    x={"value:q",axis={title="Protein half-life (mins)"}},
    y={"density:q",axis={title="Density"}},
) + quantile_plot(half_life_vals)

fig_amplitude = df |> @vlplot() + @vlplot(
    {:area, color="#a373e3"}, width=120,height=90,
    transform=[
        {density="amplitude",bandwidth=75}
    ],
    x={"value:q",axis={title="Amplitude"}},
    y={"density:q",axis={title="Density"}},
) + quantile_plot(amplitudes)

fig_period = df |> @vlplot() + @vlplot(
    {:area, color="#a373e3"}, width=120,height=90,
    transform=[
        {density="period",bandwidth=5}
    ],
    x={"value:q",axis={title="Period (mins)"}},
    y={"density:q",axis={title="Density"}},
) + quantile_plot(periods)

fig_heterogeneity = df |> fig_config + [fig_half_life fig_amplitude fig_period]
savevl(fig_heterogeneity,"output/repressilator_variability")

# Plot timeseries
(fig_deg_timeseries,df_deg_timeseries) = plot_timeseries(percentiles,sols)
savevl(fig_deg_timeseries,"output/repressilator_variability_timeseries")

## Simulate changes in ATP concentration
ATP_results_path = "output/repressilator_energy_distribution.csv"
Random.seed!(5688308)

# Load in ATP concentration data
dict_ATP_histogram = JSON.parsefile("data/Nakatani2022_ATP_histogram.json")
data_points = dict_ATP_histogram["datasetColl"][1]["data"]
_x = [d["value"][1] for d in data_points]
_y = [d["value"][2] for d in data_points]
n_points = length(_x) ÷ 2
atp_concs = [mean(_x[2*i-1:2*i]) for i in 1:n_points]
datp = (atp_concs[end]-atp_concs[1])/(n_points-1)
counts = [mean(_y[2*i-1:2*i]) for i in 1:n_points]
counts_normalised = counts ./ sum(counts)
counts_density = counts_normalised ./ datp

# Fit a normal distribution to the data
pdf_normal(x,μ,σ) = 1/(σ*sqrt(2*pi)) * exp(-(x-μ)^2/(2*σ^2))
function error_function(x)
    μ = x[1]
    σ = x[2]
    return sum((pdf_normal.(atp_concs,μ,σ) .- counts_density).^2)
end
x0 = [5.0,2.0]
res = Optim.optimize(error_function,x0,Optim.BFGS())
(μ,σ) = Optim.minimizer(res)

# Sample from truncated distribution
A_nom = 5800000
cA_nom = 9.63 # unit mM
conversion_factor = A_nom/cA_nom

d0 = Normal(μ,σ)
A_min = 10
cA_min = A_min/conversion_factor
d = truncated(d0,lower=cA_min)
n_samples = 10000
cA_vec = rand(d, n_samples)

# Run dynamic simulations on parameter distributions
@parameters A₊e
if load_results == false
    df = DataFrame(cA=Float64[], amplitude=Float64[], period=Float64[])
    @variables G1₊mRNA₊q(t) G2₊mRNA₊q(t) G3₊mRNA₊q(t) P1₊q(t) P2₊q(t) P3₊q(t) 
     for cA in cA_vec
        A = cA * conversion_factor
        ps = [A₊e => 20 + log(A/A_nom)]
        prob = ODEProblem(sys, [], tspan, ps; options...)
        sol = solve(prob, CVODE_BDF())

        (amplitude,period) = analyse_oscillations(sol)
        df_new = DataFrame(cA=cA,amplitude=amplitude,period=period)
        append!(df,df_new)
    end
    CSV.write(ATP_results_path, df)
else
    df = DataFrame(CSV.File(ATP_results_path))
end

# Extract parameters in 5, 25, 50, 75, 95th percentiles and run dynamic simuations
percentiles = [5,25,50,75,95]
cA_vals = [quantile(df.cA,p/100) for p in percentiles]
μA_vec = @. 20 + log(cA_vals*conversion_factor/A_nom)
function simulate_energy(sys,μA)
    ps = [A₊e => μA]
    prob = ODEProblem(sys, [], tspan, ps; options...)
    sol = solve(prob, CVODE_BDF())
    return sol
end

sols = [simulate_energy(sys,μr) for μr in μA_vec];
array_oscillations = analyse_oscillations.(sols)
amplitudes = [a[1] for a in array_oscillations]
periods = [a[2] for a in array_oscillations]

# Plot distributions
fig_energy_A = df |> fig_config + @vlplot(
    {:area, color="#f462a3"}, width=120,height=90,
    transform=[
        {density="cA",bandwidth=0.3}
    ],
    x={"value:q",axis={title="[ATP] (mM)"}},
    y={"density:q",axis={title="Density"}},
) + quantile_plot(cA_vals)

fig_energy_amplitide = df |> fig_config + @vlplot(
    {:area, color="#f462a3", clip=true}, width=120,height=90,
    transform=[
        {density="amplitude",bandwidth=5,minsteps=200,maxsteps=500}
    ],
    x={"value:q",axis={title="Amplitude"},scale={domain=[2100,2400]}},
    y={"density:q",axis={title="Density"}},
) + quantile_plot(amplitudes)

fig_energy_period = df |> fig_config + @vlplot(
    {:area, color="#f462a3", clip=true}, width=120,height=90,
    transform=[
        {density="period",bandwidth=0.5,minsteps=200,maxsteps=500}
    ],
    x={"value:q",axis={title="Period (mins)"},scale={domain=[50,100]}},
    y={"density:q",axis={title="Density"}},
) + quantile_plot(periods)

fig_energy = df |> fig_config + [
    fig_energy_A fig_energy_amplitide fig_energy_period
]
fig_energy.params["config"]["view"]["stroke"] = "transparent"
savevl(fig_energy,"output/repressilator_energy_distribution")

# Plot timeseries
(fig_A_timeseries,df_A_timeseries) = plot_timeseries(percentiles,sols)
savevl(fig_A_timeseries,"output/repressilator_energy_distribution_timeseries")

## Simulate variabilities in ribosome number
ribosome_results_path = "output/repressilator_ribosome_distribution.csv"
Random.seed!(5688308)

# Define lognormal distribution
cv = 0.4 # From http://book.bionumbers.org/how-much-cell-to-cell-variability-exists-in-protein-expression/
σ = sqrt(log(cv^2+1))
R_nom = 5000
d = LogNormal(log(R_nom),σ)

n_samples = 10000
R_vec = rand(d, n_samples)

# Run dynamic simulations using parameter distributions
if load_results == false
    df = DataFrame(R=Float64[], amplitude=Float64[], period=Float64[])
    @variables G1₊mRNA₊q(t) G2₊mRNA₊q(t) G3₊mRNA₊q(t) P1₊q(t) P2₊q(t) P3₊q(t) R₊q(t)
     for R in R_vec
        u0 = [R₊q => R]
        prob = ODEProblem(sys, u0, tspan, []; options...)
        sol = solve(prob, CVODE_BDF())

        (amplitude,period) = analyse_oscillations(sol)
        df_new = DataFrame(R=R,amplitude=amplitude,period=period)
        append!(df,df_new)
    end
    CSV.write(ribosome_results_path, df)
else
    df = DataFrame(CSV.File(ribosome_results_path))
end

# Extract parameters in 5, 25, 50, 75, 95th percentiles and run dynamic simuations
percentiles = [5,25,50,75,95]
R_vals = [quantile(df.R,p/100) for p in percentiles]
function simulate_ribosomes(sys,R)
    u0 = [R₊q => R]
    prob = ODEProblem(sys, u0, tspan, []; options...)
    sol = solve(prob, CVODE_BDF())
    return sol
end

sols = [simulate_ribosomes(sys,R) for R in R_vals];
array_oscillations = analyse_oscillations.(sols)
amplitudes = [a[1] for a in array_oscillations]
periods = [a[2] for a in array_oscillations]

# Plot distributions
fig_ribosome_R = df |> fig_config + @vlplot(
    {:area, color="#1aa5d3"}, width=120,height=90,
    transform=[
        {density="R",bandwidth=400}
    ],
    x={"value:q",axis={title="Ribosome number"}},
    y={"density:q",axis={title="Density"}},
) + quantile_plot(R_vals)

fig_ribosome_amplitude = df |> fig_config + @vlplot(
    {:area, color="#1aa5d3"}, width=120,height=90,
    transform=[
        {density="amplitude",bandwidth=150}
    ],
    x={"value:q",axis={title="Amplitude"}},
    y={"density:q",axis={title="Density"}},
) + quantile_plot(amplitudes)

fig_ribosome_period = df |> fig_config + @vlplot(
    {:area, color="#1aa5d3"}, width=120,height=90,
    transform=[
        {density="period",bandwidth=2}
    ],
    x={"value:q",axis={title="Period (mins)"}},
    y={"density:q",axis={title="Density"}},
) + quantile_plot(periods)

fig_ribosome = df |> fig_config + [
    fig_ribosome_R fig_ribosome_amplitude fig_ribosome_period
]
savevl(fig_ribosome,"output/repressilator_ribosome_distribution")

# Plot timeseries
(fig_R_timeseries,df_R_timeseries) = plot_timeseries(percentiles,sols)
savevl(fig_R_timeseries,"output/repressilator_ribosome_distribution_timeseries")