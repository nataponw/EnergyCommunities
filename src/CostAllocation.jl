"""
(Pending update) Analyze the contribution of respective peers via the Shapley Value method
"""
function analyzeContribution(sPeer, sY, sTS, sTec, pSca, pY, pTec, pYTec, pYTS, dT, selSolver)
    nPeer = length(sPeer)
    # Create a mapping between (coalition) id and sExcludedPeer
    map_id_sExcludedPeer = Dict{Int, Vector{Peer}}((2^nPeer-1) => Peer[])
    map_id_selection = Dict{Int, Vector{Int}}((2^nPeer-1) => ones(nPeer))
    for id ∈ 0:(2^nPeer-2)
        selection = reverse(digits(id, base=2, pad=nPeer))
        excludedPeers = sPeer[findall(x -> x == 0, selection)]
        map_id_selection[id] = selection
        map_id_sExcludedPeer[id] = excludedPeers
    end
    # Calculate the total cost of different coalition
    map_id_objvalue = Dict{Int, Float64}()
    for id ∈ 0:(2^nPeer-1)
        if sum(map_id_selection[id]) == 1
            map_id_objvalue[id] = map_id_objvalue[0]
        end
        m = initializeModel(
            sPeer, sY, sTS, sTec, pSca, pY, pTec, pYTec, pYTS, dT,
            bNoExpand=true, bIntTrade=true, sExcludedPeer=map_id_sExcludedPeer[id]
        )
        JuMP.set_optimizer(m, selSolver)
        JuMP.optimize!(m)
        map_id_objvalue[id] = JuMP.objective_value(m)
    end
    # Derive the values to the coalition
    map_id_coalvalue = Dict(0 => 0.0)
    [map_id_coalvalue[id] = (map_id_objvalue[0] - map_id_objvalue[id]) for id ∈ 1:(2^nPeer-1)];
    converter = reverse([2^x for x ∈ 0:(nPeer-1)])
    # Calculate the contribution according to the SV analysis
    map_peer_contribution = Dict{Peer, Float64}()
    for iPeer ∈ sPeer
        contribution = 0.0
        for id ∈ 0:(2^nPeer - 2)
            # Skip this coalition if iPeer ∈ selection
            if iPeer ∉ map_id_sExcludedPeer[id]
                continue
            end
            nS = sum(map_id_selection[id])
            selection_with_iPeer = map_id_selection[id] .+ Int.([x == iPeer for x ∈ sPeer])
            id_with_iPeer = sum(selection_with_iPeer .* converter)
            contribution += factorial(nS)*factorial(nPeer - nS - 1)/factorial(nPeer) * (map_id_coalvalue[id_with_iPeer] - map_id_coalvalue[id])
        end
        map_peer_contribution[iPeer] = contribution
    end
    return map_peer_contribution
end

"""
(Pending update) Analyze the cost allocation based on individual contributions and on the market equilibrium prices
"""
function analyzeCostallocation(m_cooperate, m_separate, contribution, dT)
    # Extract the indexes
    sPeer = m_cooperate[:cXC_inttrade].axes[1]
    sY = m_cooperate[:cXC_inttrade].axes[2]
    nPeer = length(sPeer); nY = length(sY)
    # Analyze the cooperative cost allocation
    costComponents = [:cXC_inttrade, :cXC_exttrade, :cCAPEX, :cDec_scap, :cElas]
    costs_separate = DataFrames.DataFrame(:peer => sPeer)
    costs_cooperate = DataFrames.DataFrame(:peer => sPeer)
    for cc ∈ costComponents
        DataFrames.insertcols!(costs_separate, cc => sum(JuMP.value.(m_separate[cc]).data, dims=2)[:])
        DataFrames.insertcols!(costs_cooperate, cc => sum(JuMP.value.(m_cooperate[cc]).data, dims=2)[:])
    end
    DataFrames.select!(costs_separate, :peer, DataFrames.AsTable(DataFrames.Not(:peer)) => sum => :tc_separate)
    DataFrames.select!(costs_cooperate, :peer, DataFrames.AsTable(DataFrames.Not(:peer)) => sum => :tc_cooperate)
    costalloc = DataFrames.outerjoin(costs_separate, costs_cooperate, on = :peer)
    costalloc = DataFrames.outerjoin(
            costalloc,
            DataFrames.DataFrame(:peer => sPeer, :contribution => [contribution[iPeer] for iPeer ∈ sPeer]),
        on = :peer)
    costalloc.tc_fairdist = costalloc.tc_separate - costalloc.contribution
    costalloc.payment_fairdist = costalloc.tc_fairdist - costalloc.tc_cooperate
    # Analze the market-based cost allocation
    mkprice_el = JuMP.dual.(m_cooperate[:ecBalac_market_el]).data
    inttrade_elImp = JuMP.value.(m_cooperate[:vXCac_inttrade_elImp]).data
    inttrade_elExp = JuMP.value.(m_cooperate[:vXCac_inttrade_elExp]).data
    mkprice_th = JuMP.dual.(m_cooperate[:ecBalac_market_th]).data
    inttrade_thImp = JuMP.value.(m_cooperate[:vXCac_inttrade_thImp]).data
    inttrade_thExp = JuMP.value.(m_cooperate[:vXCac_inttrade_thExp]).data
    payment_mkTrans_el = zeros(nPeer, nY)
    payment_mkTrans_th = zeros(nPeer, nY)
    for iPeer ∈ 1:nPeer, iY ∈ 1:nY
        payment_mkTrans_el[iPeer, iY] += dT * sum(mkprice_el[iY, :] .* inttrade_elImp[iPeer, iY, :])
        payment_mkTrans_el[iPeer, iY] -= dT * sum(mkprice_el[iY, :] .* inttrade_elExp[iPeer, iY, :])
        payment_mkTrans_th[iPeer, iY] += dT * sum(mkprice_th[iY, :] .* inttrade_thImp[iPeer, iY, :])
        payment_mkTrans_th[iPeer, iY] -= dT * sum(mkprice_th[iY, :] .* inttrade_thExp[iPeer, iY, :])
    end
    DataFrames.insertcols!(costalloc, :payment_elmk => sum(payment_mkTrans_el, dims=2)[:])
    DataFrames.insertcols!(costalloc, :payment_thmk => sum(payment_mkTrans_th, dims=2)[:])
    DataFrames.transform!(costalloc, DataFrames.AsTable([:tc_cooperate, :payment_elmk, :payment_thmk]) => sum => :tc_mkbased)
    return costalloc
end