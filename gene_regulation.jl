using BondGraphs, ModelingToolkit

## Define symbolic variables
@parameters rd α β n t r R_bS R_bP R_bI e μ0 μr μ_bI
@variables E(t)[1:3] F(t)[1:3] v(t) q(t)
D = Differential(t)

## Define new bond graph components
se = Dict(
    :description => """
    An effort source component that only uses a constant effort.

    Parameters:
        e = effort
    """,
    :numports => 1,
    :variables => Dict(
        :parameters => Dict(
            e => 0
        )
    ),
    :equations => [
        E[1] ~ e
    ]
)

re_el = Dict(
    :description => """
    A simplified model of translation elongation.

    Parameters:
        r = rate of elongation (corresponds to the reaction rate of the first translocation step)
        α = thermodynamic constant ratio between subsequent complexes = K_C_i/K_C_(i-1)
        n = number of amino acids in protein

    Ports:
        1: Ribosomal complex C_0
        2: Energetic molecule A
        3: Terminal complex C_n
    """,
    :numports => 3,
    :variables => Dict(
        :parameters => Dict(
            r => 0.0,
            α => 1,
            n => 10
        )
    ),
    :equations => [
        0 ~ v - r*(exp(E[1]) - exp(E[3] - n*E[2]))*(1-exp(E[2])/α)/((α/E[2])^n - 1),
        0 ~ F[1] - v,
        0 ~ F[2] - n*v,
        0 ~ F[3] + v
    ]
)

re_inh = Dict(
    :description => """
    A model of an enzyme with non-competitive inhibition.

    Parameters:
        r = rate of reaction
        R_bS = binding parameter for substrate
        R_bP = binding parameter for product
        R_bI = binding parameter for inhibitor

    Ports:
        1: Substrate
        2: Product
        3: Inhibitor
    """,
    :numports => 3,
    :variables => Dict(
        :parameters => Dict(
            r => 1,
            R_bS => 1,
            R_bP => 1,
            R_bI => 1
        )
    ),
    :equations => [
        0 ~ v - r*(exp(E[1]) - exp(E[2]))/(1 + exp(E[1])/R_bS + exp(E[2])/R_bP)/(1+exp(E[3])/R_bI),
        0 ~ F[1] - v,
        0 ~ F[2] + v,
        0 ~ F[3]
    ]
)

re_inh_log = Dict(
    :description => """
    A model of an enzyme with non-competitive inhibition, using a log-transformed binding parameter.

    Parameters:
        r = rate of reaction
        R_bS = binding parameter for substrate
        R_bP = binding parameter for product
        μ_bI = binding parameter for inhibitor

    Ports:
        1: Substrate
        2: Product
        3: Inhibitor
    """,
    :numports => 3,
    :variables => Dict(
        :parameters => Dict(
            r => 1,
            R_bS => 1,
            R_bP => 1,
            μ_bI => 0.0
        )
    ),
    :equations => [
        0 ~ v - r*(exp(E[1]) - exp(E[2]))/(1 + exp(E[1])/R_bS + exp(E[2])/R_bP)/(1+exp(E[3]-μ_bI)),
        0 ~ F[1] - v,
        0 ~ F[2] + v,
        0 ~ F[3]
    ]
)

ce_log = Dict(
    :description => """
    A normalised ce component with log-transformed parameter

    Parameters:
        μ0 = standard chemical potential
    """,
    :numports => 1,
    :variables => Dict(
        :parameters => Dict(
            μ0 => 0.0
        ),
        :states => Dict(
            q => 1.0
        ),
    ),
    :equations => [
        0 ~ μ0 + log(q) - E[1]
        D(q) ~ F[1]
    ]
)

re_log = Dict(
    :description => """
    A normalised ce component with log-transformed parameter

    Parameters:
        μr = reaction activation energy
    """,
    :numports => 2,
    :variables => Dict(
        :parameters => Dict(
            μr => 0.0
        )
    ),
    :equations => [
        0 ~ F[1] + F[2],
        0 ~ F[1] - (exp(μr + E[1]) - exp(μr + E[2]))
    ]
)

addlibrary!(Dict(
    :se => se, 
    :re_el => re_el, 
    :re_inh => re_inh, 
    :re_inh_log => re_inh_log, 
    :ce_log => ce_log, 
    :re_log => re_log
))

## Useful functions
"""Creates a component for a species. An SS component is used if the system is open, and a Ce component otherwise."""
function Species(name;open=true,log=false)
    if open
        return SourceSensor(;name=name)
    else
        if log
            return Component(:ce_log,name)
        else
            return Component(:ce,name)
        end
    end
end

"""Creates a component for a reaction."""
Reaction(name;log=false) = log ? Component(:re_log,name) : Component(:re,name)

modelornode(m,open) = open ? BondGraphNode(m) : m

"""Returns the parameter for a species component. The value can be either a thermodynamic constant (sym=:K) or standard potential (sym=:μ0)."""
function get_species_param(c,sym)
    if type(c) == "ce"
        if sym == :K
            return c.K
        elseif sym == :μ0
            return log(c.K)
        else
            error("Invalid symbol")
        end
    elseif type(c) == "ce_log"
        if sym == :μ0
            return c.μ0
        elseif sym == :K
            return exp(c.μ0)
        else
            error("Invalid symbol")
        end
    end
    return c
end
"""Sets the parameter for a species component. The value can be either a thermodynamic constant (sym=:K) or standard potential (sym=:μ0)."""
function set_species_param!(c,sym,val)
    if type(c) == "ce"
        if sym == :K
            c.K = val
        elseif sym == :μ0
            c.K = exp(val)
        else
            error("Invalid symbol")
        end
    elseif type(c) == "ce_log"
        if sym == :μ0
            c.μ0 = val
        elseif sym == :K
            c.μ0 = log(val)
        else
            error("Invalid symbol")
        end
    end
    return c
end
"""Sets the parameter for a species component. The value can be either a rate parameter (sym=:r) or reaction affinity (sym=:μr)."""
function set_reaction_param!(c,sym,val)
    if type(c) == "re"
        if sym == :r
            c.r = val
        elseif sym == :μr
            c.r = exp(val)
        else
            error("Invalid symbol")
        end
    elseif type(c) == "re_log"
        if sym == :μr
            c.μr = val
        elseif sym == :r
            c.μr = log(val)
        else
            error("Invalid symbol")
        end
    end
    return c
end

"""Returns a bond graph model of the enzyme catalysed reaction."""
function ECR(name="ECR";open=true)
    model = BondGraph(name)
    S = Species(:S,open=open)
    P = Species(:P,open=open)
    E = Species(:E,open=open)
    C = Component(:ce,:C)
    r1 = Component(:re,:r1)
    r2 = Component(:re,:r2)

    μE = EqualEffort(name=:μ_E)
    μC = EqualEffort(name=:μ_C)
    ES = EqualFlow(name=:ES)
    EP = EqualFlow(name=:EP)

    add_node!(model,[S,P,E,C,r1,r2,μE,μC,ES,EP])
    connect!(model,μE,E)
    connect!(model,μC,C)

    connect!(model,S,ES)
    connect!(model,μE,ES)
    connect!(model,ES,(r1,1))
    connect!(model,(r1,2),μC)

    connect!(model,μC,(r2,1))
    connect!(model,(r2,2),EP)
    connect!(model,EP,μE)
    connect!(model,EP,P)

    return modelornode(model,open)
end

"""Returns a bond graph model of transcription."""
function Transcription(name="Transcription";open=true)
    model = BondGraph(name)
    A = Species(:A; open=open)
    mRNA = Species(:mRNA; open=open)
    E = Component(:ce, :E)
    r = ECR(:ECR)
    add_node!(model,[A,mRNA,r,E])

    connect!(model,A,(r,1))
    connect!(model,(r,2),mRNA)
    connect!(model,(r,3),E)

    return modelornode(model,open)
end

"""Returns a bond graph model of transcription with an inhibitor component."""
function TranscriptionRegulated(name="Transcription";open=true)
    model = BondGraph(name)
    A = Species(:A, open=open)
    mRNA = Species(:mRNA, open=open)
    I = Species(:I, open=open)
    r = Component(:re_inh, :r)
    add_node!(model,[A,mRNA,r,I])

    connect!(model,A,(r,1))
    connect!(model,(r,2),mRNA)
    connect!(model,(r,3),I)

    return modelornode(model,open)
end

"""Returns a bond graph model of translation. The complexity of the model depends on the n_lumps parameter, which determines whether the full (n_lumps=0), simple (n_lumps=1) or multiple lumped complexes (n_lumps>1) models are used."""
function Translation(model_name="Translation";n=3,open=true,n_lumps=0,log=false)
    model = BondGraph(model_name)

    C0 = Species(:C0,open=false)
    C0_potential = EqualEffort(name=:μ_C0)
    Cx_name = (n_lumps>0) ? :B0 : :C1
    Cx = Component(:ce,Cx_name)
    Cx_potential = EqualEffort(name=Symbol("μ_$Cx_name"))
    Cn_name = (n_lumps>0) ? :Cn : Symbol("C$n")
    Cn = Component(:ce_log,Cn_name)
    Cn_potential = EqualEffort(name=Symbol("μ_$Cn_name"))

    M = Species(:M,open=open)
    M_potential = EqualEffort(name=:μ_M)
    R = Species(:R,open=open)
    R_potential = EqualEffort(name=:μ_R)
    A = Species(:A,open=open)
    A_potential = EqualEffort(name=:μ_A)
    P = Species(:P,open=open)

    MR = EqualFlow(name=:MR)
    MRP = EqualFlow(name=:MRP)

    association = Component(:re,:rb)
    r1 = Component(:re,:r1)
    r1f = EqualFlow(name=:r1f)
    r1r = EqualFlow(name=:r1r)
    if n_lumps == 0
        r = Elongation(:elong,m=n-1,log=log)
    elseif n_lumps == 1
        r = Component(:re_el,:elong)
        r.n = n-1
    else
        r = ElongationBlocks(:elong,m=n-1,n_lumps=n_lumps)
    end
    termination = Component(:re_log,:rd)

    add_node!(model,[
        Cx;
        Cx_potential;
        Cn;
        Cn_potential;
        M;
        M_potential;
        R;
        R_potential;
        A;
        A_potential;
        P;
        MR;
        MRP;
        association;
        r;
        termination;
        C0;
        C0_potential;
        r1;
        r1f;
        r1r
    ])

    connect!(model,C0_potential,C0)
    connect!(model,Cx_potential,Cx)
    connect!(model,Cn_potential,Cn)
    connect!(model,A_potential,A)
    connect!(model,M_potential,M)
    connect!(model,R_potential,R)

    # Elongation
    connect!(model,Cx_potential,(r,1))
    connect!(model,A_potential,(r,2))
    connect!(model,(r,3),Cn_potential)

    # Termination
    connect!(model,Cn_potential,(termination,1))
    connect!(model,(termination,2),MRP)
    #connect!(model,MRP,M_potential)
    connect!(model,MRP,R_potential)
    connect!(model,MRP,P)

    # Ribosome complexation
    connect!(model,M_potential,MR)
    connect!(model,R_potential,MR)
    connect!(model,MR,(association,1))
    connect!(model,(association,2),C0_potential)

    # First elongation step
    connect!(model,C0_potential,r1f)
    connect!(model,A_potential,r1f)
    connect!(model,r1f,(r1,1))
    connect!(model,(r1,2),r1r)
    connect!(model,r1r,Cx_potential)
    connect!(model,r1r,M_potential)

    return modelornode(model,open)
end

"""Returns the full bond graph model of elongation."""
function Elongation(model_name="Elongation";m=3,open=true,log=false)
    model = BondGraph(model_name)

    X0 = Species(:X0,open=open,log=log)
    Xm = Species("X$m",open=open,log=log)
    intermediates = [Species("X$i",open=false,log=log) for i in 1:m-1]
    complexes = [intermediates;Xm]
    complex_potentials = [EqualEffort(name="μ_$(name(c))") for c in complexes]

    A = Species(:A,open=open,log=log)
    A_potential = EqualEffort(name=:μ_A)

    reactions = [Reaction("r$i";log=log) for i in 1:m]
    reaction_junctions = [EqualFlow(name="v_$(name(r))") for r in reactions]

    add_node!(model,[
        X0;
        A;
        A_potential
        complexes;
        complex_potentials;
        reactions;
        reaction_junctions;
    ])

    # Connect species to junctions
    for (c,μ) in zip(complexes,complex_potentials)
        connect!(model,μ,c)
    end
    connect!(model,A_potential,A)

    # Elongation
    for (i,r) in enumerate(reactions)
        if i == 1
            connect!(model,X0,reaction_junctions[i])
        else
            connect!(model,complex_potentials[i-1],reaction_junctions[i])
        end
        connect!(model,A_potential,reaction_junctions[i])
        connect!(model,reaction_junctions[i],(r,1))
        connect!(model,(r,2),complex_potentials[i])
    end

    return modelornode(model,open)
end

"""Returns the number of complexes required in each lump for a model with multiple lumped complexes."""
function lumped_complex_lengths(n,n_lumps)
    l_lumps = floor(Int,n/n_lumps)
    r = n%n_lumps
    m_lumps = ones(Int,n_lumps)*l_lumps
    for i in 1:r
        m_lumps[i] += 1
    end
    return m_lumps
end

"""Returns a bond graph model of elongation with multiple lumped complexes."""
function ElongationBlocks(model_name="Elongation";m=19,n_lumps=5,open=true)
    n_internal_lumps = n_lumps-1
    model = BondGraph(model_name)

    X0 = Species(:X0,open=open)
    Xm = Species("X$m",open=open)
    intermediates = [Species("B$i",open=false) for i in 1:n_internal_lumps]
    complexes = [intermediates;Xm]
    complex_potentials = [EqualEffort(name="μ_$(name(c))") for c in complexes]

    A = Species(:A,open=open)
    A_potential = EqualEffort(name=:μ_A)

    reactions = [Component(:re_el,"el$i") for i in 1:n_lumps]

    m_lumps = lumped_complex_lengths(m,n_lumps)
    for (i,r) in enumerate(reactions)
        r.n = m_lumps[i]
    end

    add_node!(model,[
        X0;
        A;
        A_potential
        complexes;
        complex_potentials;
        reactions;
    ])

    # Connect species to junctions
    for (c,μ) in zip(complexes,complex_potentials)
        connect!(model,μ,c)
    end
    connect!(model,A_potential,A)

    # Elongation
    for (i,r) in enumerate(reactions)
        if i == 1
            connect!(model,X0,(r,1))
        else
            connect!(model,complex_potentials[i-1],(r,1))
        end
        connect!(model,A_potential,(r,2))
        connect!(model,(r,3),complex_potentials[i])
    end

    return modelornode(model,open)
end

"""Returns a bond graph model of degradation."""
function Degradation(name="deg";open=true)
    model = BondGraph(name)

    S = Species(:S,open=open)
    deg = Component(:re_log,:deg)
    P = Component(:se,:P)
    add_node!(model,[S,deg,P])

    connect!(model,S,(deg,1))
    connect!(model,(deg,2),P)

    return modelornode(model,open)
end

"""Returns a bond graph model of gene expression."""
function GeneExpression(name="GeneExpression"; n=1200, open=true, log=true, n_lumps=1)
    model = BondGraph(name)

    A = Species(:A; open=open)
    I = Species(:I, open=open)
    R = Species(:R; open=open)
    P = Species(:P; open=open,log=true)
    mRNA = Component(:ce,:mRNA)

    μmRNA = EqualEffort(name=:μmRNA)
    μA = EqualEffort(name=:μA)
    μP = EqualEffort(name=:μP)

    Tc = Component(:re_inh_log,:Tc)
    Tl = Translation(:Tl,n=n;n_lumps=n_lumps,log=log)
    mRNA_d = Degradation(:mRNA_deg)
    P_d = Degradation(:P_deg)

    add_node!(model,[A,I,R,P,mRNA,μmRNA,μA,Tc,Tl,mRNA_d,μP,P_d])

    # Connect species to potential junctions
    connect!(model,μmRNA,mRNA)
    connect!(model,μA,A)
    connect!(model,μP,P)

    # Transcription
    connect!(model,μA,(Tc,1))
    connect!(model,(Tc,2),μmRNA)
    connect!(model,I,(Tc,3))

    # Translation
    connect!(model,μmRNA,(Tl,1))
    connect!(model,(Tl,2),R)
    connect!(model,μA,(Tl,3))
    connect!(model,(Tl,4),μP)

    # Degradation
    connect!(model,μmRNA,mRNA_d)
    connect!(model,μP,P_d)

    return modelornode(model,open)
end