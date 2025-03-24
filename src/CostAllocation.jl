"""
    iterate_over_all_coalition(m)

Given a model `m`, solve the problem over all possible coalitions and collect the costs
"""
function iterate_over_all_coalition(m)
    # Prepare the list of coalition setup
    sPeer = m[:sPeer]; sY = m[:sY]; sTS = m[:sTS]
    nPeer = length(sPeer)
    dc_idx_sExcludedPeer = Dict{Int, Vector{Peer}}()
    dc_idx_selection = Dict{Int, Vector{Int}}()
    for id ∈ 0:(2^nPeer-1)
        selection = digits(id, base=2, pad=nPeer)
        excludedPeers = sPeer[findall(x -> x == 0, selection)]
        dc_idx_selection[id] = selection
        dc_idx_sExcludedPeer[id] = excludedPeers
    end
    # Iterate through each coalition
    dc_idx_objvalue = Dict{Int, Float64}()
    dc_idx_costs = Dict{Int, NamedTuple}()
    for id ∈ 0:(2^nPeer-1)
        if sum(dc_idx_selection[id]) == 1
            dc_idx_objvalue[id] = dc_idx_objvalue[0]
            dc_idx_costs[id] = dc_idx_costs[0]
            continue
        end
        # Reverse existing peer exclusion
        for iPeer ∈ sPeer
            if JuMP.is_fixed(m[:vXCac_inttrade_elImp][iPeer, first(sY), first(sTS)])
                JuMP.unfix.(m[:vXCac_inttrade_elImp][iPeer, sY, sTS])
                JuMP.unfix.(m[:vXCac_inttrade_elExp][iPeer, sY, sTS])
                JuMP.unfix.(m[:vXCac_inttrade_thImp][iPeer, sY, sTS])
                JuMP.unfix.(m[:vXCac_inttrade_thExp][iPeer, sY, sTS])
                JuMP.set_lower_bound.(m[:vXCac_inttrade_elImp][iPeer, sY, sTS], 0.0)
                JuMP.set_lower_bound.(m[:vXCac_inttrade_elExp][iPeer, sY, sTS], 0.0)
                JuMP.set_lower_bound.(m[:vXCac_inttrade_thImp][iPeer, sY, sTS], 0.0)
                JuMP.set_lower_bound.(m[:vXCac_inttrade_thExp][iPeer, sY, sTS], 0.0)
            end
        end
        # Apply new peer exclusion
        JuMP.fix.(m[:vXCac_inttrade_elImp][dc_idx_sExcludedPeer[id], sY, sTS], 0.0, force=true)
        JuMP.fix.(m[:vXCac_inttrade_elExp][dc_idx_sExcludedPeer[id], sY, sTS], 0.0, force=true)
        JuMP.fix.(m[:vXCac_inttrade_thImp][dc_idx_sExcludedPeer[id], sY, sTS], 0.0, force=true)
        JuMP.fix.(m[:vXCac_inttrade_thExp][dc_idx_sExcludedPeer[id], sY, sTS], 0.0, force=true)
        # Optimize the model with new peer exclusion
        JuMP.optimize!(m)
        dc_idx_objvalue[id] = JuMP.objective_value(m)
        dc_idx_costs[id] = (
            cXC_inttrade = JuMP.value.(m[:cXC_inttrade]),
            cXC_exttrade = JuMP.value.(m[:cXC_exttrade]),
            cCAPEX = JuMP.value.(m[:cCAPEX]),
            cOPEX = JuMP.value.(m[:cOPEX]),
            cDec_scap = JuMP.value.(m[:cDec_scap]),
            cElas = JuMP.value.(m[:cElas]),
            cFS = JuMP.value.(m[:cFS]),
        )
    end
    return (; dc_idx_objvalue, dc_idx_costs)
end

"""
    allocate_cost_contribution(dc_idx_objvalue, dc_idx_costs, sPeer)

Given the output of `iterate_over_all_coalition`, calculate the expected (final) total cost (`tc_final`) and corresponding internal payment (`chipin`) based on individual contributions
"""
function allocate_cost_contribution(dc_idx_objvalue, dc_idx_costs, sPeer)
    # Calculate the contribution
    benefits = [dc_idx_objvalue[id] for id ∈ 0:(2^length(sPeer)-1)]
    benefits = first(benefits) .- benefits
    contribution = Dict(zip(sPeer, analyze_shapleyvalue(benefits)))
    # Calculate the ex post `chipin`
    allocatedCost = DataFrames.DataFrame(
        peer = String[],
        tc_ext_nocoop = Float64[],
        tc_external = Float64[],
        contribution = Float64[]
    )
    idx_nocoop = 0; idx_wicoop = 2^length(sPeer)-1
    for p ∈ sPeer
        tc_ext_nocoop = 0.; tc_external = 0.
        for costobj ∈ keys(dc_idx_costs[idx_nocoop])
            tc_ext_nocoop += sum(dc_idx_costs[idx_nocoop][costobj][p, :])
            tc_external += sum(dc_idx_costs[idx_wicoop][costobj][p, :])
        end
        push!(allocatedCost, (p.value, tc_ext_nocoop, tc_external, contribution[p]))
    end
    DataFrames.transform!(allocatedCost, [:tc_ext_nocoop, :contribution] => (-) => :tc_final)
    DataFrames.transform!(allocatedCost, [:tc_final, :tc_external] => (-) => :chipin)
    return allocatedCost
end

"""
    analyze_shapleyvalue(benefits::Vector{Float64})

Given `benefits`, the list of benefits from cooperation, calculate individual contribution.

The location indicates the coalition setup, see `calculate_contribution_chipin`. The first element is with an empty coalition, and the last element is with the grand coalition.
"""
function analyze_shapleyvalue(benefits::Vector{Float64})
    nP = Int(log2(length(benefits)))
    map_idx_participation = [digits(Bool, idx, base=2, pad=nP) for idx ∈ 0:(length(benefits)-1)]
    converter = [2^i for i ∈ 0:(nP-1)]
    contribution = zeros(nP)
    for iP ∈ 1:nP, idx ∈ 1:2^nP
        if !map_idx_participation[idx][iP]
            nS = sum(map_idx_participation[idx])
            idx_with_iP = idx + converter[iP]
            contribution[iP] += factorial(nS)*factorial(nP - nS - 1)/factorial(nP) * (benefits[idx_with_iP] - benefits[idx])
        end
    end
    return contribution
end

"""
    allocate_cost_marginalprice(m::JuMP.Model)

Given an optimized model `m`, process the cost allocation based on the marginal prices
"""
function allocate_cost_marginalprice(m::JuMP.Model)
    # Extract market prices
    price_mk_el = JuMP.dual.(m[:ecBalac_market_el])
    price_mk_th = JuMP.dual.(m[:ecBalac_market_th])
    # Extract cost components
    allocatedCost = DataFrames.DataFrame(peer=String[], year=Int[], variable=String[], value=Float64[])
    for p ∈ m[:sPeer], y ∈ m[:sY]
        ## External costs
        push!(allocatedCost, (p.value, y.value, "cXC_exttrade", JuMP.value.(m[:cXC_exttrade])[p, y]))
        push!(allocatedCost, (p.value, y.value, "cXC_inttrade", JuMP.value.(m[:cXC_inttrade])[p, y]))
        push!(allocatedCost, (p.value, y.value, "cCAPEX", JuMP.value.(m[:cCAPEX])[p, y]))
        push!(allocatedCost, (p.value, y.value, "cOPEX", JuMP.value.(m[:cOPEX])[p, y]))
        push!(allocatedCost, (p.value, y.value, "cDec_scap", JuMP.value.(m[:cDec_scap])[p, y]))
        push!(allocatedCost, (p.value, y.value, "cElas", JuMP.value.(m[:cElas])[p, y]))
        push!(allocatedCost, (p.value, y.value, "cFS", JuMP.value.(m[:cFS])[p, y]))
        ## Internal energy trade, payment based on the marginal prices and fees
        push!(allocatedCost, (p.value, y.value, "inttrade_el_mk", m[:dT] * (
            + sum(price_mk_el[y, :] .* JuMP.value.(m[:vXCac_inttrade_elImp])[p, y, :])
            - sum(price_mk_el[y, :] .* JuMP.value.(m[:vXCac_inttrade_elExp])[p, y, :])
        )))
        push!(allocatedCost, (p.value, y.value, "inttrade_th_mk", m[:dT] * (
            + sum(price_mk_th[y, :] .* JuMP.value.(m[:vXCac_inttrade_thImp])[p, y, :])
            - sum(price_mk_th[y, :] .* JuMP.value.(m[:vXCac_inttrade_thExp])[p, y, :])
        )))
    end
    # Process into the fial structure
    allocatedCost = DataFrames.combine(DataFrames.groupby(allocatedCost, [:peer, :variable]), :value => sum => :value)
    allocatedCost = DataFrames.unstack(allocatedCost, :peer, :variable, :value)
    DataFrames.select!(allocatedCost, :peer, [:cXC_exttrade, :cXC_inttrade, :cCAPEX, :cOPEX, :cDec_scap, :cElas, :cFS] => (+) => :tc_external, :inttrade_el_mk, :inttrade_th_mk)
    DataFrames.transform!(allocatedCost, :, [:tc_external, :inttrade_el_mk, :inttrade_th_mk] => (+) => :tc_final)
    return allocatedCost
end
