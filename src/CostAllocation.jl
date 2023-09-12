"""
    iterate_over_all_coalition(m)

Given a model `m`, solve the problem over all possible coalition
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
        for iPeer ∈ sPeer, iY ∈ sY
            if JuMP.is_fixed(m[:vXCac_inttrade_elImp][iPeer, iY, first(sTS)])
                JuMP.unfix.(m[:vXCac_inttrade_elImp][iPeer, iY, sTS])
                JuMP.unfix.(m[:vXCac_inttrade_elExp][iPeer, iY, sTS])
                JuMP.unfix.(m[:vXCac_inttrade_thImp][iPeer, iY, sTS])
                JuMP.unfix.(m[:vXCac_inttrade_thExp][iPeer, iY, sTS])
                JuMP.set_lower_bound.(m[:vXCac_inttrade_elImp][iPeer, iY, sTS], 0.0)
                JuMP.set_lower_bound.(m[:vXCac_inttrade_elExp][iPeer, iY, sTS], 0.0)
                JuMP.set_lower_bound.(m[:vXCac_inttrade_thImp][iPeer, iY, sTS], 0.0)
                JuMP.set_lower_bound.(m[:vXCac_inttrade_thExp][iPeer, iY, sTS], 0.0)
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
            cDec_scap = JuMP.value.(m[:cDec_scap]),
            cElas = JuMP.value.(m[:cElas]),
            cFS = JuMP.value.(m[:cFS]),
        )
    end
    return (; dc_idx_objvalue, dc_idx_costs)
end

"""
    calculate_chipin_coorperative(dc_idx_objvalue, dc_idx_costs, sPeer)

Given the output of `iterate_over_all_coalition`, calculate the chipin based on coorperative contributions
"""
function calculate_contribution_chipin(dc_idx_objvalue, dc_idx_costs, sPeer)
    # Calculate the contribution
    benefits = [dc_idx_objvalue[id] for id ∈ 0:(2^length(sPeer)-1)]
    benefits = first(benefits) .- benefits
    contribution = Dict(zip(sPeer, ESAAnalytics.shapleyvalueanalysis(benefits)))
    # Calculate the ex post `chipin`
    costs = DataFrames.DataFrame(
        peer = String[],
        tc_ext_nocoop = Float64[],
        tc_ext_wicoop = Float64[],
        contribution = Float64[]
    )
    idx_nocoop = 0; idx_wicoop = 2^length(sPeer)-1
    for p ∈ sPeer
        tc_ext_nocoop = 0.; tc_ext_wicoop = 0.
        for costobj ∈ keys(dc_idx_costs[idx_nocoop])
            tc_ext_nocoop += sum(dc_idx_costs[idx_nocoop][costobj][p, :])
            tc_ext_wicoop += sum(dc_idx_costs[idx_wicoop][costobj][p, :])
        end
        push!(costs, (p.value, tc_ext_nocoop, tc_ext_wicoop, contribution[p]))
    end
    DataFrames.transform!(costs, [:tc_ext_nocoop, :contribution] => (-) => :tc_expected)
    DataFrames.transform!(costs, [:tc_expected, :tc_ext_wicoop] => (-) => :chipin)
    return costs
end
