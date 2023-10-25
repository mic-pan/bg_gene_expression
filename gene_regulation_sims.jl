# Option to run full models (which could take a long time)
run_full_model = false

using DifferentialEquations, Sundials, Plots, DataFrames, VegaLite, CSV
include("gene_circuits.jl")
include("vl_functions.jl")
include("chemical_systems.jl")

## Run simulation of simple gene expression model
# Simulation options for simple model
n_lumps = 1
simple = (n_lumps==1)
options = Dict(
    :solver => CVODE_BDF(),
    :abstol => 1e-12,
    :reltol => 1e-9
)

# Define model of gene expression
model = GeneExpression(open=false,n_lumps=n_lumps)
A = Component(:se,:A)
swap!(model,model.A,A)
set_params_gr!(model,simple=simple)
I = Component(:Se,:I)
swap!(model,model.I,I)
μ0_P = model.P.μ0
f(t) = (t<0.0) ? μ0_P + log(5000) : μ0_P + log(1e-3)
@register_symbolic f(t)
model.I.es = (t -> f(t))

# Run simulation
tspan = (-250.0,500.0)
sol = simulate(model,tspan;options...);

# Extract results into DataFrame
@variables mRNA₊q(t) P₊q(t)
t_plot = -100:1:400
df_gr = DataFrame(
    t=t_plot,
    I=model.I.es.(t_plot) .- model.Tc.μ_bI,
    M=sol(t_plot,idxs=mRNA₊q).u,
    P=sol(t_plot,idxs=P₊q).u,
    model=:GR,
    gene=1
)

# Plot results
fig_gr_I = df_gr |> vl_config() + @vlplot(
    {:line, color="#1657a8"}, width=120, height=90,
    x={:t, axis={title=""}}, 
    y={:I, axis={title="Inhibition potential"}}
)
fig_gr_M = df_gr |> vl_config() + @vlplot(
    {:line, color="#1657a8"}, width=120, height=90,
    x={:t, axis={title=""}}, 
    y={:M, axis={title="mRNA"}}
)
fig_gr_P = df_gr |> vl_config() + @vlplot(
    {:line, color="#1657a8"}, width=120, height=90,
    x={:t, axis={title="Time (min)"}}, 
    y={:P, axis={title="Protein"}}
)
fig_gr = df_gr |> vl_config() + [fig_gr_I;fig_gr_M;fig_gr_P]
fig_gr.params["config"]["view"]["stroke"] = :transparent

# Plot energetics of the model
# Transcription
@variables Tc₊p1₊E(t) Tc₊p2₊E(t) Tc₊p3₊E(t) Tc₊p1₊F(t) Tc₊p2₊F(t) Tc₊p3₊F(t)
P_Tc = Tc₊p1₊E*Tc₊p1₊F + Tc₊p2₊E*Tc₊p2₊F + Tc₊p3₊E*Tc₊p3₊F
P_Tc_vec = sol(t_plot,idxs=P_Tc).u
df_Tc = DataFrame(t=t_plot,P=P_Tc_vec,process="Transcription")

# Complex binding
@variables Tl₊rb₊p1₊E(t) Tl₊rb₊p2₊E(t) Tl₊rb₊p1₊F(t) Tl₊rb₊p2₊F(t)
P_b = Tl₊rb₊p1₊E*Tl₊rb₊p1₊F + Tl₊rb₊p2₊E*Tl₊rb₊p2₊F
P_b_vec = sol(t_plot,idxs=P_b).u
df_b = DataFrame(t=t_plot,P=P_b_vec,process="Complex binding")

# First elongation reaction
@variables Tl₊r1₊p1₊E(t) Tl₊r1₊p2₊E(t) Tl₊r1₊p1₊F(t) Tl₊r1₊p2₊F(t)
P_1 = Tl₊r1₊p1₊E*Tl₊r1₊p1₊F + Tl₊r1₊p2₊E*Tl₊r1₊p2₊F

# Translation elongation
@variables Tl₊elong₊p1₊E(t) Tl₊elong₊p2₊E(t) Tl₊elong₊p3₊E(t) Tl₊elong₊p1₊F(t) Tl₊elong₊p2₊F(t) Tl₊elong₊p3₊F(t)
P_elong = Tl₊elong₊p1₊E*Tl₊elong₊p1₊F + Tl₊elong₊p2₊E*Tl₊elong₊p2₊F + Tl₊elong₊p3₊E*Tl₊elong₊p3₊F
P_elong_tot = P_1 + P_elong
P_elong_tot_vec = sol(t_plot,idxs=P_elong_tot).u
df_elong = DataFrame(t=t_plot,P=P_elong_tot_vec,process="Translation elongation")

# Translation termination
@variables Tl₊rd₊p1₊E(t) Tl₊rd₊p2₊E(t) Tl₊rd₊p1₊F(t) Tl₊rd₊p2₊F(t)
P_t = Tl₊rd₊p1₊E*Tl₊rd₊p1₊F + Tl₊rd₊p2₊E*Tl₊rd₊p2₊F
P_t_vec = sol(t_plot,idxs=P_t).u
df_t = DataFrame(t=t_plot,P=P_t_vec,process="Translation termination")

# mRNA degradation
@variables mRNA_deg₊deg₊p1₊E(t) mRNA_deg₊deg₊p2₊E(t) mRNA_deg₊deg₊p1₊F(t) mRNA_deg₊deg₊p2₊F(t)
P_Mdeg = mRNA_deg₊deg₊p1₊E*mRNA_deg₊deg₊p1₊F + mRNA_deg₊deg₊p2₊E*mRNA_deg₊deg₊p2₊F
P_Mdeg_vec = sol(t_plot,idxs=P_Mdeg).u
df_Mdeg = DataFrame(t=t_plot,P=P_Mdeg_vec,process="mRNA degradation")

# Protein degradation
@variables P_deg₊deg₊p1₊E(t) P_deg₊deg₊p2₊E(t) P_deg₊deg₊p1₊F(t) P_deg₊deg₊p2₊F(t)
P_Pdeg = P_deg₊deg₊p1₊E*P_deg₊deg₊p1₊F + P_deg₊deg₊p2₊E*P_deg₊deg₊p2₊F
P_Pdeg_vec = sol(t_plot,idxs=P_Pdeg).u
df_Pdeg = DataFrame(t=t_plot,P=P_Pdeg_vec,process="Protein degradation")

df_power = [df_Tc; df_b; df_elong; df_t; df_Mdeg; df_Pdeg]

fig_power = df_power |> vl_config() + @vlplot(
    :line, width=240, height=180,
    x={:t, axis={title="Time (mins)"}}, 
    y={:P, axis={title="Power (kT/min)",grid=true},scale={type=:log}, grid=true},
    color={
        :process, 
        scale={scheme=:category10},
        sort=["Transcription","Complex binding","Translation elongation",
            "Translation termination", "mRNA degradation", "Protein degradation"]
    }
)

## Run simulations of the toggle switch
# Define toggle switch model
model = ToggleSwitch(n=1200,n_lumps=n_lumps,log=true)
set_params_toggle!(model,n=1200,simple=simple)

# Run simulations on toggle switch model
@variables G1₊mRNA₊q(t) P1₊q(t) G2₊mRNA₊q(t) P2₊q(t)
tspan = (0.0,200.0)
sys = ODESystem(model)
prob = ODEProblem(sys,[P1₊q => 50.0, P2₊q => 100.0],tspan,[];options...)
sol_toggle1 = solve(prob,CVODE_BDF())
prob = ODEProblem(sys,[P1₊q => 200.0, P2₊q => 100.0],tspan,[];options...)
sol_toggle2 = solve(prob,CVODE_BDF())

# Extract results into DataFrame
function make_df_toggle(sol,index)
    t_plot = 0:1:200
    df_G1 = DataFrame(
        t = t_plot,
        M = sol(t_plot,idxs=G1₊mRNA₊q).u,
        P = sol(t_plot,idxs=P1₊q).u,
        model = :Toggle,
        sim = index,
        gene = 1
    )
    df_G2 = DataFrame(
        t = t_plot,
        M = sol(t_plot,idxs=G2₊mRNA₊q).u,
        P = sol(t_plot,idxs=P2₊q).u,
        model = :Toggle,
        sim = index,
        gene = 2
    )
    df_toggle_sim1 = [df_G1;df_G2]
end
df1 = make_df_toggle(sol_toggle1,1)
df2 = make_df_toggle(sol_toggle2,2)
df_toggle = [df1;df2]

# Plot results
fig_toggle_M = df_toggle |> vl_config() + @vlplot(
    :line, width=120, height=90,
    x={:t, axis={title=""}}, 
    y={:M, axis={title="mRNA"}},
    color={
        :gene,
        type=:ordinal,
        scale={
            range=["#1657a8", "#048059"]
        }
    },
    strokeDash={:sim}
)

fig_toggle_P = df_toggle |> vl_config() + @vlplot(
    :line, width=120, height=90,
    x={:t, axis={title="Time (min)"}}, 
    y={:P, axis={title="Protein"}},
    color={
        :gene,
        type=:ordinal,
        scale={
            range=["#1657a8", "#048059"]
        }
    },
    strokeDash={:sim}
)

fig_toggle = df_toggle |> vl_config() + [fig_toggle_M;fig_toggle_P]
fig_toggle.params["config"]["view"]["stroke"] = :transparent

## Run simulations of repressilator model
# Construct model of repressilator
model = Repressilator(n=1200,n_lumps=n_lumps,log=true)
set_params_repressilator!(model,n=1200,simple=simple)

# Run simulations
tspan = (0.0,500.0)
sol = simulate(model,tspan;options...)

# Compile results into DataFrame
@variables G1₊mRNA₊q(t) G2₊mRNA₊q(t) G3₊mRNA₊q(t) P1₊q(t) P2₊q(t) P3₊q(t) R₊q(t)
function make_df_repressilator(sol)
    t_plot = 0.0:1.0:150.0
    df_repressilator1 = DataFrame(
        t = t_plot,
        M = sol(t_plot,idxs=G1₊mRNA₊q).u,
        P = sol(t_plot,idxs=P1₊q).u,
        model = :Repressilator,
        gene = "Gene 1"
    )
    df_repressilator2 = DataFrame(
        t = t_plot,
        M = sol(t_plot,idxs=G2₊mRNA₊q).u,
        P = sol(t_plot,idxs=P2₊q).u,
        model = :Repressilator,
        gene = "Gene 2"
    )
    df_repressilator3 = DataFrame(
        t = t_plot,
        M = sol(t_plot,idxs=G3₊mRNA₊q).u,
        P = sol(t_plot,idxs=P3₊q).u,
        model = :Repressilator,
        gene = "Gene 3"
    )
    df_repressilator = [df_repressilator1;df_repressilator2;df_repressilator3]
end
df_repressilator = make_df_repressilator(sol)

# Plot results
fig_repressilator_M = df_repressilator |> vl_config() + @vlplot(
    :line, width=120, height=90,
    x={:t, axis={title=""}}, 
    y={:M, axis={title="mRNA"}},
    color={
        :gene,
        type=:ordinal,
        scale={
            range=["#1657a8", "#048059", "#6626bb"]
        }
    }
)

fig_repressilator_P = df_repressilator |> vl_config() + @vlplot(
    :line, width=120, height=90,
    x={:t, axis={title="Time (min)"}}, 
    y={:P, axis={title="Protein"}},
    color={
        :gene,
        type=:ordinal,
        scale={
            range=["#1657a8", "#048059", "#6626bb"]
        }
    },
)

fig_repressilator = df_repressilator |> vl_config() + [fig_repressilator_M;fig_repressilator_P]
fig_repressilator.params["config"]["view"]["stroke"] = :transparent

savevl(fig_gr,"output/gr_sim")
savevl(fig_power,"output/gr_sim_power")
savevl(fig_toggle,"output/toggle_sim")
savevl(fig_repressilator,"output/repressilator_sim")

## Compare reduced models to full models
# File names
gr_sim_data = "output/gr_sim_full.csv"
toggle_sim_data = "output/toggle_sim_full.csv"
repressilator_sim_data = "output/repressilator_sim_full.csv"

# Run simulations of full models
n_lumps = 0
simple = (n_lumps==1)

if run_full_model == true
    # Model of gene expression
    model = GeneExpression(open=false,n_lumps=n_lumps)
    A = Component(:se,:A)
    swap!(model,model.A,A)
    set_params_gr!(model,simple=simple)
    I = Component(:Se,:I)
    swap!(model,model.I,I)
    μ0_P = model.P.μ0
    f(t) = (t<0.0) ? μ0_P + log(5000) : μ0_P + log(1e-3)
    @register_symbolic f(t)
    model.I.es = (t -> f(t))

    tspan = (-250.0,500.0)
    sys = chemical_odesys(model)
    prob = ODEProblem(sys,[],tspan,[];options...)
    sol = solve(prob,CVODE_BDF())

    t_plot = -100:1:400
    df_gr_full = DataFrame(
        t=t_plot,
        I=model.I.es.(t_plot) .- model.Tc.μ_bI,
        M=sol(t_plot,idxs=mRNA₊q).u,
        P=sol(t_plot,idxs=P₊q).u,
        model=:GR,
        gene=1
    )
    CSV.write(gr_sim_data,df_gr_full)

    # Model of toggle switch
    model = ToggleSwitch(n=1200,n_lumps=n_lumps,log=true)
    set_params_toggle!(model,n=1200,simple=simple)
    tspan = (0.0,200.0)
    sys = chemical_odesys(model)
    prob = ODEProblem(sys,[P1₊q => 50.0, P2₊q => 100.0],tspan,[];options...)
    sol_toggle1 = solve(prob,CVODE_BDF())
    prob = ODEProblem(sys,[P1₊q => 200.0, P2₊q => 100.0],tspan,[];options...)
    sol_toggle2 = solve(prob,CVODE_BDF())
    
    df1 = make_df_toggle(sol_toggle1,1)
    df2 = make_df_toggle(sol_toggle2,2)
    df_toggle_full = [df1;df2]
    CSV.write(toggle_sim_data,df_toggle_full)

    # Model of repressilator
    model = Repressilator(n=1200,n_lumps=n_lumps,log=true)
    set_params_repressilator!(model,n=1200,simple=simple)
    tspan = (0.0,500.0)
    odesys = chemical_odesys(model)
    prob = ODEProblem(odesys,[],tspan,[];options...)
    sol = solve(prob,CVODE_BDF())
    df_repressilator_full = make_df_repressilator(sol)
    CSV.write(repressilator_sim_data,df_repressilator)
end

# Plot comparison for gene regulation
df_gr_full = DataFrame(CSV.File(gr_sim_data))
df_gr2 = copy(df_gr)
df_gr2[!,:type] .= "Reduced"
df_gr_full2 = copy(df_gr_full)
df_gr_full2[!,:type] .= "Full"
df_gr_comparison = vcat(df_gr_full2,df_gr2)

fig_gr_M_c = df_gr_comparison |> vl_config() + @vlplot(
    {:line, color="#1657a8"}, width=120, height=90,
    x={:t, axis={title=""}}, 
    y={:M, axis={title="mRNA"}},
    color={:type, scale={range=["#1657a8", "#f462a3"]}},
    strokeDash={:type, scale={range=[[4,0],[4,4]]}}
)

fig_gr_P_c = df_gr_comparison |> vl_config() + @vlplot(
    {:line, color="#1657a8"}, width=120, height=90,
    x={:t, axis={title="Time (min)"}}, 
    y={:P, axis={title="Protein"}},
    color={:type, scale={range=["#1657a8", "#f462a3"]}},
    strokeDash={:type, scale={range=[[4,0],[4,4]]}}
)

fig_gr_c = df_gr_comparison |> vl_config() + [fig_gr_M_c;fig_gr_P_c]
fig_gr_c.params["config"]["view"]["stroke"] = :transparent

# Plot comparison for toggle switch
df_toggle_full = DataFrame(CSV.File(toggle_sim_data))
df_toggle2 = subset(df_toggle, :sim => ByRow(x -> x == 1))
df_toggle2.type .= "Reduced"
df_toggle_full2 = subset(df_toggle_full, :sim => ByRow(x -> x == 1))
df_toggle_full2.type .= "Full"
df_toggle_comparison = vcat(df_toggle_full2,df_toggle2) 

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

fig_toggle_M_c = fig_config + @vlplot(
    :line, data=df_toggle_full2, width=120, height=90,
    x={"t:q", axis={title=""}}, 
    y={"M:q", axis={title="mRNA"}},
    color={
        :gene,
        type=:ordinal,
        scale={
            range=["#1657a8", "#048059"]
        }
    }
) + @vlplot(
    {:line, color="#f462a3", strokeDash=[4,4]}, data=df_toggle2, width=120, height=90,
    x={"t:q", axis={title=""}}, 
    y={"M:q", axis={title="mRNA"}},
    color={
        :gene,
        type=:ordinal,
        scale={
            range=["#f462a3"]
        }
    }
)

fig_toggle_P_c = fig_config + @vlplot(
    :line, data=df_toggle_full2, width=120, height=90,
    x={"t:q", axis={title=""}}, 
    y={"P:q", axis={title="Protein"}},
    color={
        :gene,
        type=:ordinal,
        scale={
            range=["#1657a8", "#048059"]
        }
    }
) + @vlplot(
    {:line, color="#f462a3", strokeDash=[4,4]}, data=df_toggle2, width=120, height=90,
    x={"t:q", axis={title="Time (min)"}}, 
    y={"P:q", axis={title="Protein"}},
    color={
        :gene,
        type=:ordinal,
        scale={
            range=["#f462a3"]
        }
    }
)

# Plot comparison for repressilator
df_repressilator_full = DataFrame(CSV.File(repressilator_sim_data))
df_repressilator2 = copy(df_repressilator)
df_repressilator2.type .= "Reduced"
df_repressilator_full2 = copy(df_repressilator_full)
df_repressilator_full2.type .= "Full"

fig_repressilator_M_c = fig_config + @vlplot(
    :line, data=df_repressilator_full2, width=120, height=90,
    x={"t:q", axis={title=""}}, 
    y={"M:q", axis={title="mRNA"}},
    color={
        :gene,
        type=:ordinal,
        scale={
            range=["#1657a8", "#048059", "#6626bb"]
        }
    }
) + @vlplot(
    {:line, color="#f462a3", strokeDash=[4,4]}, data=df_repressilator2, width=120, height=90,
    x={"t:q", axis={title=""}}, 
    y={"M:q", axis={title="mRNA"}},
    color={
        :gene,
        type=:ordinal,
        scale={
            range=["#f462a3"]
        }
    }
)

fig_repressilator_P_c = fig_config + @vlplot(
    :line, data=df_repressilator_full2, width=120, height=90,
    x={"t:q", axis={title="Time (min)"}}, 
    y={"P:q", axis={title="Protein"}},
    color={
        :gene,
        type=:ordinal,
        scale={
            range=["#1657a8", "#048059", "#6626bb"]
        }
    }
) + @vlplot(
    {:line, color="#f462a3", strokeDash=[4,4]}, data=df_repressilator2, width=120, height=90,
    x={"t:q", axis={title=""}}, 
    y={"P:q", axis={title="Protein"}},
    color={
        :gene,
        type=:ordinal,
        scale={
            range=["#f462a3"]
        }
    }
)

# Save comparison figures
savevl(fig_gr_c,"output/gr_sim_comparison")
savevl(fig_toggle_M_c,"output/toggle_sim_comparison_M")
savevl(fig_toggle_P_c,"output/toggle_sim_comparison_P")
savevl(fig_repressilator_M_c,"output/repressilator_sim_comparison_M")
savevl(fig_repressilator_P_c,"output/repressilator_sim_comparison_P")