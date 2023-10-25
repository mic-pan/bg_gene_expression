using ModelingToolkit, BondGraphs, SparseArrays

"""Returns the acausal relationships encoded by a bond graph."""
function get_constraints(bg)
    (subsystems, connections) = BondGraphs.get_subsys_and_connections(bg)
    _sys = BondGraphs.compose_bg_model(subsystems, connections, bg.name, false)
    sys = expand_connections(_sys)
    return sys
end

extract_name(x) = string(x)[1:end-3]
isflow(x) = (x==:F || x==:v)

"""Classifies a variable defined as a string. Returns the variable type."""
classify(x::String) = return Symbol(x[end])
"""Classifies a variable defined as a symbolic. Returns the variable type."""
classify(x::Term) = classify(extract_name(x))
"""Returns a classification of the efforts, flows and stores in an ODE system."""
function classify(sys::ODESystem)
    s = states(sys)
    var_names = extract_name.(s)
    var_types = classify.(var_names)
    efforts = [x for (x,v) in zip(s,var_types) if v == :E]
    flows = [x for (x,v) in zip(s,var_types) if v == :F]
    stores = [x for (x,v) in zip(s,var_types) if v == :q]
    return (efforts,flows,stores)
end
"""Classifies an equation based on the variables it relates."""
function classify(eq::Equation)
    vars = [x for x in get_variables(eq) if x isa Term]
    var_types = classify.(vars)
    if all(var_types .== :E)
        return :effort
    elseif all(isflow.(var_types))
        return :flow
    elseif (:q in var_types) && any(isflow.(var_types))
        return :differential
    elseif (:q in var_types) && (:E in var_types)
        return :chemical_potential
    elseif (:E in var_types) && any(isflow.(var_types))
        return :reaction
    else
        return :other
    end
end

"""Determines the algebraic form of an effort"""
function solve_cp(eq)
    vars = get_dynamic_vars(eq)
    var_types = classify.(vars)
    i = findfirst(var_types .== :E)
    e = vars[i]
    sol = ModelingToolkit.solve_for(eq,e)
    return (e,sol)
end
"""Returns substitutions required for efforts"""
function cp_substitutions(cp_eqs)
    tuples = solve_cp.(cp_eqs)
    cp_efforts = [x[1] for x in tuples]
    cp_sols = [x[2] for x in tuples]
    return cp_efforts,cp_sols
end
get_dynamic_vars(eq) = [x for x in get_variables(eq) if x isa Term]

"""Determines the algebraic form of a flow"""
function solve_r(eq)
    vars = get_dynamic_vars(eq)
    var_types = classify.(vars)
    i = findfirst(isflow.(var_types))
    f = vars[i]
    sol = ModelingToolkit.solve_for(eq,f)
    return (f,sol)
end
"""Returns substitutions required for flows"""
function r_substitutions(r_eqs)
    tuples = solve_r.(r_eqs)
    cp_flows = [x[1] for x in tuples]
    cp_sols = [x[2] for x in tuples]
    return cp_flows,cp_sols
end
r_efforts(r_eqs) = [x for x in vcat(get_dynamic_vars.(r_eqs)...) if classify(x) == :E]
ismember(x,y) = any(isequal.(x,y))

"""Returns substitutions required for efforts in a system"""
function solve_reaction_potentials(e_eqs,cp_subs)
    cp_vars = cp_subs[1]
    all_efforts = unique(vcat(get_dynamic_vars.(e_eqs)...))
    unknown_vars = [x for x in all_efforts if !ismember(x,cp_vars)]
    sol = ModelingToolkit.solve_for(e_eqs,unknown_vars)
    return Dict(zip(unknown_vars,sol))
end
"""Returns substitutions required for flows in a system"""
function solve_species_flows(f_eqs,r_subs)
    r_flows = r_subs[1]
    all_flows = unique(vcat(get_dynamic_vars.(f_eqs)...))
    unknown_vars = [x for x in all_flows if !ismember(x,r_flows)]
    sol = ModelingToolkit.solve_for(f_eqs,unknown_vars)
    return Dict(zip(unknown_vars,sol))
end
"""
Returns ODEs for a biochemical bond graph. The algorithm assumes that bond graphs are only composed of species and reactions, which allows it to be faster than more general procedures.
"""
function chemical_odesys(bg)
    sys = get_constraints(bg)
    # Classify equations
    eqs = equations(sys)
    eq_types = classify.(eqs)

    cp_eqs = [eq for (eq,t) in zip(eqs,eq_types) if t == :chemical_potential]
    e_eqs = [eq for (eq,t) in zip(eqs,eq_types) if t == :effort]
    r_eqs = [eq for (eq,t) in zip(eqs,eq_types) if t == :reaction]
    f_eqs = [eq for (eq,t) in zip(eqs,eq_types) if t == :flow]
    d_eqs = [eq for (eq,t) in zip(eqs,eq_types) if t == :differential]

    # Solve for reaction potentials
    cp_subs = cp_substitutions(cp_eqs)
    effort_subs = solve_linear_sys(e_eqs,cp_subs)

    # Solve for species flows
    r_subs = r_substitutions(r_eqs)
    flow_subs = solve_linear_sys(f_eqs,r_subs)

    # Substitute variables
    d_eqs1 = substitute(d_eqs,flow_subs)
    d_eqs2 = substitute(d_eqs1,Dict(zip(r_subs...)))
    d_eqs3 = substitute(d_eqs2,effort_subs)
    d_eqs4 = substitute(d_eqs3,Dict(zip(cp_subs...)))

    @parameters t
    reduced_sys = ODESystem(d_eqs4, t; name=Symbol(name(bg)), defaults=sys.defaults)
    return reduced_sys
end

""" Solves a system of symbolic matrix equations. The function decomposes the equations into symbolic and numerical components to speed up the solver."""
function solve_linear_sys(eqs,subs)
    vars = subs[1]
    all_efforts = unique(vcat(get_dynamic_vars.(eqs)...))
    unknown_vars = [x for x in all_efforts if !ismember(x,vars)]

    @assert length(unknown_vars) == length(eqs)
    unknown_reverse_index = Dict(x => (i,:A) for (i,x) in enumerate(unknown_vars))
    input_reverse_index = Dict(x => (i,:b) for (i,x) in enumerate(vars))
    reverse_index = merge(unknown_reverse_index,input_reverse_index)

    A = spzeros(Int,length(eqs),length(unknown_vars))
    b = zeros(Int,length(eqs),length(vars))
    c = zeros(Num,length(eqs))
    
    for (i,eq) in enumerate(eqs)
        rhs = eq.rhs
        lhs = eq.lhs
        update!(A,b,c,rhs,i,1,reverse_index)
        update!(A,b,c,lhs,i,-1,reverse_index)
    end

    d,additional_terms = process_terms(c)

    sol = (A\b)*vars + (A\d)*additional_terms
    return Dict(zip(unknown_vars,sol))
end
function add_to_matrices!(A,b,c,k,v,i,sign,reverse_index)
    if haskey(reverse_index,k)
        (j,flag) = reverse_index[k]
        if flag == :A
            A[i,j] += sign*v
        elseif flag == :b
            b[i,j] -= sign*v
        end
    else
        c[i] -= sign*v*k
    end
end
function update!(A,b,c,expr::Symbolics.Add,i,sign,reverse_index)
    for (k,v) in expr.dict
        add_to_matrices!(A,b,c,k,v,i,sign,reverse_index)
    end
end
function update!(A,b,c,expr,i,sign,reverse_index)
    add_to_matrices!(A,b,c,expr,1,i,sign,reverse_index)
end
function process_terms(c)
    additional_terms = Num[]
    idx_additional = Int[]
    for (i,_x) in enumerate(c)
        if !iszero(_x)
            push!(additional_terms,_x)
            push!(idx_additional,i)
        end
    end
    d = zeros(Int,length(c),length(idx_additional))
    for (i,j) in enumerate(idx_additional)
        d[j,i] = 1
    end
    return d,additional_terms
end