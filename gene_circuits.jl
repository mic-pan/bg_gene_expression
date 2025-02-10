include("gene_regulation.jl")

"""Add specific component for transformer used in regulation"""
TF2 = Dict(
    :description => "Linear Transformer with modulus 2",
    :numports => 2,
    :equations => [
        0 ~ E[2] - 2 * E[1],
        0 ~ F[1] + 2 * F[2]
    ],
)

TF4 = Dict(
    :description => "Linear Transformer with modulus 2",
    :numports => 2,
    :equations => [
        0 ~ E[2] - 4 * E[1],
        0 ~ F[1] + 4 * F[2]
    ],
)

addlibrary!(Dict(
    :TF2 => TF2,
    :TF4 => TF4,
))

## Model parameters
n = 1200 # Corresponding to an amino acid length of 300
α = exp(1.25)
μ_folding = 20.0
μ0_P = n*log(α) - μ_folding

## General functions
"""Calculate thermodynamic constant from amount and chemical potential."""
compute_K(x,μ) = (1/x)*exp(μ)

"""Set parameters for transcription model"""
function set_params_Tc!(model,w,θ,K_A,μ0_I)
    model.r = w/(K_A*θ)
    model.R_bS = K_A*θ
    model.R_bP = 1e6
    model.μ_bI = μ0_I + log(30) # Achieves inhibition when I > 30
    return nothing
end

"""Set parameters for simple translation model"""
function set_params_Tl!(model,n,α,kf,μA_nom)
    model.C0.K = 1.0
    model.C0.q = 1e-6

    model.r1.r = kf/exp(μA_nom)
    K_C1 = α
    model.B0.K = K_C1/(n-1)
    model.B0.q = 1e-6*(n-1)

    model.Cn.μ0 = n*log(α)
    model.Cn.q = 1e-6

    model.rb.r = 1e-3
    model.elong.r = kf/exp(μA_nom)/K_C1*α
    model.elong.n = n-1
    model.elong.α = α
    model.rd.μr = log(100*kf) - model.Cn.μ0
    return nothing
end

"""Set parameters for full translation model"""
function set_params_Tl_full!(model,n,α,kf,μA_nom)
    set_species_param!(model.C0,:μ0,0.0)
    model.C0.q = 1e-6

    set_species_param!(model.C1,:μ0,log(α))
    model.C1.q = 1e-6
    model.r1.r = kf/exp(μA_nom)
    set_species_param!(model.r1,:μr,log(kf)-μA_nom)
    intermediates = Dict(
        i => getproperty(model.elong,Symbol("X$(i-1)")) for i in 2:n-1
    )
    intermediates[1] = model.C1
    intermediates[n] = getproperty(model,Symbol("C$n"))

    for i in 2:n
        μ0_prev = get_species_param(intermediates[i-1],:μ0)
        c = intermediates[i]
        set_species_param!(c,:μ0,μ0_prev + log(α))
        c.q = 1e-6
        r = getproperty(model.elong,Symbol("r$(i-1)"))
        set_reaction_param!(r, :μr, log(kf) - μ0_prev - μA_nom)
    end

    set_reaction_param!(model.rd, :μr, log(100*kf) - get_species_param(intermediates[n],:μ0))
    model.rb.r = 1e-3

    return model
end

""" Set parameters for gene expression model"""
function set_params_gr!(model;n=n, A_nom=5800000, μA_nom=20, simple=true)
    K_A = compute_K(A_nom,μA_nom)
    #K_P = α^n/log(μ_folding)

    if model isa BondGraph
        model.A.e = μA_nom
        model.R.K = 1.0
        model.R.q = 5000.0 
        model.P.μ0 = μ0_P
        model.P.q = 1e-6
    end

    model.mRNA.K = 1.0
    model.mRNA.q = 1e-6

    w = 4.14
    θ = 4.38
    set_params_Tc!(model.Tc,w,θ,K_A,μ0_P)
    
    γ_max = 1260
    kf = 4*γ_max # Multiplied by 4 to account for 4 ATP molecules per amino acid
    if simple
        set_params_Tl!(model.Tl,n,α,kf,μA_nom)
    else
        set_params_Tl_full!(model.Tl,n,α,kf,μA_nom)
    end

    model.mRNA_deg.deg.μr  = log(0.1) # Todo: consider changing to ln(2)/2, consistent with synthetic circuits of Weisse et al.
    model.mRNA_deg.P.e = log(1e-6) # Small value corresponding to nucleotide potential
    model.P_deg.deg.μr = log(0.02) - μ0_P # Todo: consider changing to ln(2)/4, consistent with synthetic circuits of Weisse et al.
    model.P_deg.P.e = 0.0 # Zero, corresponding to approximate potential of amino acid monomers
    return nothing
end

## Define the model of the toggle switch
"""Returns a bond graph model of the toggle switch."""
function ToggleSwitch(name="ToggleSwitch";n=1200,n_lumps=1,log=true,h=2)
    model = BondGraph(name)
    A = Component(:se,:A)
    μA = EqualEffort(name=:μA)
    R = Component(:ce,:R)
    μR = EqualEffort(name=:μR)

    G1 = GeneExpression(:G1,n=n,n_lumps=n_lumps,log=log)
    P1 = Component(:ce_log,:P1)
    μP1 = EqualEffort(name=:μP1)

    G2 = GeneExpression(:G2,n=n,n_lumps=n_lumps,log=log)
    P2 = Component(:ce_log,:P2)
    μP2 = EqualEffort(name=:μP2)

    tf_symb = Symbol("TF$h")
    TF1 = Component(tf_symb,:h1)
    TF2 = Component(tf_symb,:h2)

    add_node!(model,[A,μA,R,μR,G1,P1,μP1,G2,P2,μP2,TF1,TF2])

    # Common potentials
    connect!(model,μA,A)
    connect!(model,μR,R)
    connect!(model,μP1,P1)
    connect!(model,μP2,P2)

    # Shared resources
    connect!(model,(G1,1),μA)
    connect!(model,μR,(G1,3))
    connect!(model,(G2,1),μA)
    connect!(model,μR,(G2,3))

    # Protein synthesis
    connect!(model,(G1,4),μP1)
    connect!(model,(G2,4),μP2)
    
    # Inhibition
    connect!(model,μP1,(TF1,1))
    connect!(model,μP2,(TF2,1))
    connect!(model,(TF2,2),(G1,2))
    connect!(model,(TF1,2),(G2,2))
    
    return model
end

"""Set parameters for the toggle switch."""
function set_params_toggle!(model;n=1200,simple=true,h=2)
    set_params_gr!(model.G1,n=n,simple=simple)
    set_params_gr!(model.G2,n=n,simple=simple)

    model.A.e = 20.0
    model.R.K = 1.0
    model.R.q = 5000.0

    model.G1.Tl.rb.r = 1e-2
    model.G2.Tl.rb.r = 1e-2

    μ0_P = n*log(α) - μ_folding
    model.P1.μ0 = μ0_P
    model.P1.q = 1e-6
    model.P2.μ0 = μ0_P
    model.P2.q = 1e-6

    model.G1.Tc.μ_bI = h*(μ0_P + log(100)) 
    model.G2.Tc.μ_bI = h*(μ0_P + log(100))

    model.G1.mRNA_deg.deg.μr = log(log(2)/2)
    model.G2.mRNA_deg.deg.μr = log(log(2)/2)

    model.G1.P_deg.deg.μr = log(log(2)/4) - μ0_P
    model.G2.P_deg.deg.μr = log(log(2)/4) - μ0_P
    return nothing
end

## Define the model of the repressilator
"""Returns a bond graph model of the repressilator"""
function Repressilator(name="Repressilator";n=1200,n_lumps=1,log=true,h=2)
    model = BondGraph(name)
    A = Component(:se,:A)
    μA = EqualEffort(name=:μA)
    R = Component(:ce,:R)
    μR = EqualEffort(name=:μR)

    G1 = GeneExpression(:G1,n=n,n_lumps=n_lumps,log=log)
    P1 = Component(:ce_log,:P1)
    μP1 = EqualEffort(name=:μP1)

    G2 = GeneExpression(:G2,n=n,n_lumps=n_lumps,log=log)
    P2 = Component(:ce_log,:P2)
    μP2 = EqualEffort(name=:μP2)

    G3 = GeneExpression(:G3,n=n,n_lumps=n_lumps,log=log)
    P3 = Component(:ce_log,:P3)
    μP3 = EqualEffort(name=:μP3)

    tf_symb = Symbol("TF$h")
    TF3 = Component(tf_symb,:h3)
    TF1 = Component(tf_symb,:h1)
    TF2 = Component(tf_symb,:h2)

    add_node!(model,[A,μA,R,μR,G1,P1,μP1,G2,P2,μP2,G3,P3,μP3,TF1,TF2,TF3])

    # Common potentials
    connect!(model,μA,A)
    connect!(model,μR,R)
    connect!(model,μP1,P1)
    connect!(model,μP2,P2)
    connect!(model,μP3,P3)

    # Shared resources
    connect!(model,(G1,1),μA)
    connect!(model,μR,(G1,3))
    connect!(model,(G2,1),μA)
    connect!(model,μR,(G2,3))
    connect!(model,(G3,1),μA)
    connect!(model,μR,(G3,3))

    # Protein synthesis
    connect!(model,(G1,4),μP1)
    connect!(model,(G2,4),μP2)
    connect!(model,(G3,4),μP3)
    
    # Inhibition
    connect!(model,μP1,(TF1,1))
    connect!(model,μP2,(TF2,1))
    connect!(model,μP3,(TF3,1))
    connect!(model,(TF2,2),(G1,2))
    connect!(model,(TF3,2),(G2,2))
    connect!(model,(TF1,2),(G3,2))
    
    return model
end

"""Set parameters for the repressilator model."""
function set_params_repressilator!(model;n=1200,simple=true,h=2)
    A_nom = 5800000
    μA = 20

    set_params_gr!(model.G1,n=n,simple=simple)
    set_params_gr!(model.G2,n=n,simple=simple)
    set_params_gr!(model.G3,n=n,simple=simple)

    model.A.e = μA
    model.R.K = 1.0
    model.R.q = 5000.0

    model.G1.Tl.rb.r = 1e-2
    model.G2.Tl.rb.r = 1e-2
    model.G3.Tl.rb.r = 1e-2

    μ0_P = n*log(α) - μ_folding
    model.P1.μ0 = μ0_P
    model.P1.q = 100.0
    model.P2.μ0 = μ0_P
    model.P2.q = 100.0
    model.P3.μ0 = μ0_P
    model.P3.q = 1000.0

    model.G1.Tc.μ_bI = h*(μ0_P + log(100))
    model.G2.Tc.μ_bI = h*(μ0_P + log(100))
    model.G3.Tc.μ_bI = h*(μ0_P + log(100))

    model.G1.mRNA_deg.deg.μr = log(log(2)/2)
    model.G2.mRNA_deg.deg.μr = log(log(2)/2)
    model.G3.mRNA_deg.deg.μr = log(log(2)/2)

    model.G1.P_deg.deg.μr = log(log(2)/4) - μ0_P
    model.G2.P_deg.deg.μr = log(log(2)/4) - μ0_P
    model.G3.P_deg.deg.μr = log(log(2)/4) - μ0_P

    return nothing
end