"""
    initialize_model(sPeer, sY, sTS, pSca, pY, pTec, pYTec, pYTS, dT; kwargs)

Build an optimization problem

# Arguments
- `sPeer`, `sY`, `sTS`, `sTec` : sets of peers, years, timesteps, and technologies, respectively.
- `pSca`, `pY`, `pTec`, `pYTec`, `pYTS` : Extracted parameter dataframes from [`process_parameters`](@ref)
- `dT` : Time interval

# Keyword Arguments
- `bOneoff`: for each peer, each technology can be installed once over the investigated years.
- `bFuelswitch`: the community can either import natural gas or hydrogen, but not both.
- `bConElas`: activate the elastic thermal demand module. Note that in doing so, the model becomes (MI)QP.
- `bNoExpand`: deactivate the expansion and decommission modules, i.e., only the starting capacities are considered.
- `bCHPThCurtail`: allow the curtailment of CHPs' thermal generation
- `bIntTrade`: activate the internal trading module, otherwise, the community is considered as a single entity.
- `sExcludedPeer`: List of peers excluded from the internal trading
- `solverbackend`: if a solver is provided, the model is created via `JuMP.direct_model` instead of `JuMP.Model`
"""
function initialize_model(
    sPeer, sY, sTS, sTec, pSca, pY, pTec, pYTec, pYTS, dT;
    bOneoff::Bool=false, bFuelswitch::Bool=false, bConElas::Bool=false, bNoExpand::Bool=false,
    bCHPThCurtail::Bool=false, bIntTrade::Bool=false, sExcludedPeer=Peer[],
    solverbackend=missing,
    )
    # Sorting indexes and parameters
    sort!(sPeer); sort!(sY); sort!(sTS); sort!(sTec)
    sort!(pSca, [:agentid])
    sort!(pY, [:year, :agentid])
    sort!(pTec, [:tec, :agentid])
    sort!(pYTec, [:tec, :year, :agentid])
    sort!(pYTS, [:timestep, :year, :agentid])
    # Building the model
    if ismissing(solverbackend)
        m = JuMP.Model()
    else
        m = JuMP.direct_model(solverbackend)
    end
    # Essential mappings and parameters
    sTS_next = Dict(sTS .=> circshift(sTS, -1))
    timesteps_in_year = Int(hours_in_year/dT)
    nPeer = length(sPeer); nY = length(sY); nTS = length(sTS); nTec = length(sTec)
    # Precalculation --------------------------------------------------------------
    # Calculate annuity rates for new investment
    DataFrames.transform!(pTec, [:ncap_loaninterest, :ncap_loanduration] => (
        (r, n) -> r ./ (1 .- (1 .+ r).^(-n))
    ) => :ncap_annuityrate)
    # Attach essential constants and sets
    m[:dT] = dT
    m[:sPeer] = sPeer; m[:sY] = sY; m[:sTS] = sTS; m[:sTec] = sTec
    # Attach essential parameters
    m[:dem_el_hh] = DenseAxisArray(reshape(pYTS.dem_el_hh, nPeer, nY, nTS), sPeer, sY, sTS)
    m[:dem_el_ev] = DenseAxisArray(reshape(pYTS.dem_el_ev, nPeer, nY, nTS), sPeer, sY, sTS)
    m[:dem_th_hh] = DenseAxisArray(reshape(pYTS.dem_th_hh, nPeer, nY, nTS), sPeer, sY, sTS)
    m[:impprice_el] = DenseAxisArray(reshape(pY.impprice_el, nPeer, nY), sPeer, sY)
    m[:expprice_el] = DenseAxisArray(reshape(pY.expprice_el, nPeer, nY), sPeer, sY)
    m[:impprice_th] = DenseAxisArray(reshape(pY.impprice_th, nPeer, nY), sPeer, sY)
    m[:expprice_th] = DenseAxisArray(reshape(pY.expprice_th, nPeer, nY), sPeer, sY)
    m[:impprice_ng] = DenseAxisArray(reshape(pY.impprice_ng, nPeer, nY), sPeer, sY)
    m[:impprice_h2] = DenseAxisArray(reshape(pY.impprice_h2, nPeer, nY), sPeer, sY)
    m[:inttrade_fee_el] = DenseAxisArray(reshape(filter(:parameter => ==("inttrade_fee_el"), pSca).value, nPeer), sPeer)
    m[:inttrade_fee_th] = DenseAxisArray(reshape(filter(:parameter => ==("inttrade_fee_th"), pSca).value, nPeer), sPeer)
    # Energy exchange variables ---------------------------------------------------
    ## Agent
    JuMP.@variable(m, 0.0 ≤ vXCph_peer_elImp[sPeer, sY, sTS])
    JuMP.@variable(m, 0.0 ≤ vXCph_peer_elExp[sPeer, sY, sTS])
    JuMP.@variable(m, 0.0 ≤ vXCph_peer_thImp[sPeer, sY, sTS])
    JuMP.@variable(m, 0.0 ≤ vXCph_peer_thExp[sPeer, sY, sTS])
    ## Community
    JuMP.@variable(m, 0.0 ≤ vXCph_comm_elImp[sY, sTS])
    JuMP.@variable(m, 0.0 ≤ vXCph_comm_elExp[sY, sTS])
    JuMP.@variable(m, 0.0 ≤ vXCph_comm_thImp[sY, sTS])
    JuMP.@variable(m, 0.0 ≤ vXCph_comm_thExp[sY, sTS])
    JuMP.@variable(m, 0.0 ≤ vXCph_comm_ngImp[sY, sTS])
    JuMP.@variable(m, 0.0 ≤ vXCph_comm_h2Imp[sY, sTS])
    # Investment and capacity variables -------------------------------------------
    JuMP.@variable(m, 0.0 ≤ vCap[iPeer ∈ sPeer, sY, iTec ∈ sTec] ≤ gp(pTec, :pot, [iPeer.value, iTec.value]))
    JuMP.@variable(m, 0.0 ≤ vCap_new[sPeer, sY, sTec])
    JuMP.@variable(m, 0.0 ≤ vInv[sPeer, sY, sTec])
    JuMP.@variable(m, zDec_scap[sPeer, sY, sTec], Bin)
    # Operation variables ---------------------------------------------------------
    ## Electricity
    JuMP.@variable(m, 0.0 ≤ vOpt_pv_gen[sPeer, sY, sTS])
    JuMP.@variable(m, 0.0 ≤ vOpt_wt_gen[sPeer, sY, sTS])
    JuMP.@variable(m, 0.0 ≤ vOpt_es_ch[sPeer, sY, sTS])
    JuMP.@variable(m, 0.0 ≤ vOpt_es_dch[sPeer, sY, sTS])
    JuMP.@variable(m, 0.0 ≤ vOpt_es_soc[sPeer, sY, sTS])
    JuMP.@variable(m, 0.0 ≤ vOpt_chpng_gen_el[sPeer, sY, sTS])
    JuMP.@variable(m, 0.0 ≤ vOpt_chph2_gen_el[sPeer, sY, sTS])
    ## Heat
    JuMP.@variable(m, 0.0 ≤ vOpt_gb_gen[sPeer, sY, sTS])
    JuMP.@variable(m, 0.0 ≤ vOpt_hb_gen[sPeer, sY, sTS])
    JuMP.@variable(m, 0.0 ≤ vOpt_hp_gen[sPeer, sY, sTS])
    JuMP.@variable(m, 0.0 ≤ vOpt_ts_ch[sPeer, sY, sTS])
    JuMP.@variable(m, 0.0 ≤ vOpt_ts_dch[sPeer, sY, sTS])
    JuMP.@variable(m, 0.0 ≤ vOpt_ts_soc[sPeer, sY, sTS])
    JuMP.@variable(m, 0.0 ≤ vOpt_chpng_gen_th[sPeer, sY, sTS])
    JuMP.@variable(m, 0.0 ≤ vOpt_chph2_gen_th[sPeer, sY, sTS])
    ## Fuel consumption
    JuMP.@variable(m, 0.0 ≤ vCon_gb_ng[sPeer, sY, sTS])
    JuMP.@variable(m, 0.0 ≤ vCon_hb_h2[sPeer, sY, sTS])
    JuMP.@variable(m, 0.0 ≤ vCon_hp_el[sPeer, sY, sTS])
    JuMP.@variable(m, 0.0 ≤ vCon_chpng_ng[sPeer, sY, sTS])
    JuMP.@variable(m, 0.0 ≤ vCon_chph2_h2[sPeer, sY, sTS])
    # Demand elasticity variables -------------------------------------------------
    ## Demand shift : dem_el_ev
    JuMP.@variable(m, vDem_evshift_increase[sPeer, sY, sTS])
    JuMP.@variable(m, vDem_evshift_soc[sPeer, sY, sTS])
    ## Thermal curtailment
    JuMP.@variable(m, 0.0 ≤ vDem_thCur[sPeer, sY, sTS])
    JuMP.@variable(m, 0.0 ≤ vDem_thCur_rate[sPeer, sY])
    # Troubleshooting and failsafe variables --------------------------------------
    JuMP.@variable(m, 0.0 ≤ vFS_sink_el[sPeer, sY, sTS])
    JuMP.@variable(m, 0.0 ≤ vFS_source_el[sPeer, sY, sTS])
    JuMP.@variable(m, 0.0 ≤ vFS_sink_th[sPeer, sY, sTS])
    JuMP.@variable(m, 0.0 ≤ vFS_source_th[sPeer, sY, sTS])
    # Capacity constraints --------------------------------------------------------
    ## New capacity and investment costs
    JuMP.@constraint(m, ecInv,
        vInv .== reshape(pYTec.inv, nPeer, nY, nTec) .* vCap_new
    )
    ## Installed capacity logic
    par_scap_cap = DenseAxisArray(reshape(pTec.scap_cap, nPeer, nTec), sPeer, sTec)
    par_scap_finalyear = DenseAxisArray(reshape(pTec.scap_finalyear, nPeer, nTec), sPeer, sTec)
    par_ncap_lifetime = DenseAxisArray(reshape(pTec.ncap_lifetime, nPeer, nTec), sPeer, sTec)
    JuMP.@constraint(m, ecCap[iPeer ∈ sPeer, iY ∈ sY, iTec ∈ sTec],
        + vCap[iPeer, iY, iTec]
        ==
        + sum(vCap_new[iPeer, iiY, iTec] for iiY ∈ sY if (iiY ≤ iY) && (iY < iiY + par_ncap_lifetime[iPeer, iTec]))
        + par_scap_cap[iPeer, iTec] *
        (iY ≤ par_scap_finalyear[iPeer, iTec] ? 1.0 : 0.0) * # Enforce lifetime of starting capacities
        (1.0 - sum(zDec_scap[iPeer, iiY, iTec] for iiY ∈ sY if iiY ≤ iY)) # Early decommission of starting capacities
    )
    ## Prevent multiple decommision which can act as early decommission of new capacities
    JuMP.@constraint(m, icDec_scap[iPeer ∈ sPeer, iTec ∈ sTec],
        sum(zDec_scap[iPeer, :, iTec]) ≤ 1.0
    )
    # Operation constraints -------------------------------------------------------
    ## PV: generation and curtailment
    JuMP.@constraint(m, icOpt_pv,
        vOpt_pv_gen .≤ reshape(pYTS.gen_pv, nPeer, nY, nTS) .* repeat(vCap[:, :, Tec("pv")], inner=(1, 1, nTS))
    )
    ## WT: generation and curtailment
    JuMP.@constraint(m, icOpt_wt,
        vOpt_wt_gen .≤ reshape(pYTS.gen_wt, nPeer, nY, nTS) .* repeat(vCap[:, :, Tec("wt")], inner=(1, 1, nTS))
    )
    ## ES: maximum charge, discharge, and SOC; operation (SOC-loop)
    JuMP.@constraint(m, icOpt_es1,
        vOpt_es_ch .+ vOpt_es_dch .≤ repeat(vCap[:, :, Tec("es")], inner=(1, 1, nTS))
    )
    JuMP.@constraint(m, icOpt_es2,
        vOpt_es_soc .≤ repeat(vCap[:, :, Tec("es")], inner=(1, 1, nTS))
    )
    par_eff_st_es = DenseAxisArray(filter(:tec => ==("es"), pTec).eff_st, sPeer)
    par_eff_ch_es = DenseAxisArray(filter(:tec => ==("es"), pTec).eff_ch, sPeer)    
    JuMP.@constraint(m, ecOpt_es[iPeer ∈ sPeer, iY ∈ sY, iTS ∈ sTS],
        + vOpt_es_soc[iPeer, iY, sTS_next[iTS]]
        ==
        + vOpt_es_soc[iPeer, iY, iTS] * par_eff_st_es[iPeer]^(dT/hours_in_day)
        + dT*(
            + vOpt_es_ch[iPeer, iY, iTS] * sqrt(par_eff_ch_es[iPeer])
            - vOpt_es_dch[iPeer, iY, iTS] / sqrt(par_eff_ch_es[iPeer])
        )
    )
    ## GB: thermal generation and fuel consumption
    JuMP.@constraint(m, icOpt_gb,
        vOpt_gb_gen .≤ repeat(vCap[:, :, Tec("gb")], inner=(1, 1, nTS))
    )
    JuMP.@constraint(m, ecCon_gb,
        vOpt_gb_gen .== vCon_gb_ng .* repeat(filter(:tec => ==("gb"), pTec).eff_th, inner=(1, nY, nTS))
    )
    ## HB: thermal generation and fuel consumption
    JuMP.@constraint(m, icOpt_hb,
        vOpt_hb_gen .≤ repeat(vCap[:, :, Tec("hb")], inner=(1, 1, nTS))
    )
    JuMP.@constraint(m, ecCon_hb,
        vOpt_hb_gen .== vCon_hb_h2 .* repeat(filter(:tec => ==("hb"), pTec).eff_th, inner=(1, nY, nTS))
    )
    ## HP: thermal generation and electricity consumption
    JuMP.@constraint(m, icOpt_hp,
        vOpt_hp_gen .≤ repeat(vCap[:, :, Tec("hp")], inner=(1, 1, nTS))
    )
    JuMP.@constraint(m, ecCon_hp,
        vOpt_hp_gen .== vCon_hp_el .* reshape(pYTS.cop, nPeer, nY, nTS)
    )
    ## DH: thermal import and export
    JuMP.@constraint(m, icOpt_dh,
        vXCph_peer_thImp .+ vXCph_peer_thExp .≤ repeat(vCap[:, :, Tec("dh")], inner=(1, 1, nTS))
    )
    ## TS: maximum charge, discharge, and SOC; operation (SOC-loop)
    JuMP.@constraint(m, icOpt_ts1,
        vOpt_ts_ch .+ vOpt_ts_dch .≤ repeat(vCap[:, :, Tec("ts")], inner=(1, 1, nTS))
    )
    JuMP.@constraint(m, icOpt_ts2,
        vOpt_ts_soc .≤ repeat(vCap[:, :, Tec("ts")], inner=(1, 1, nTS))
    )
    par_eff_st_ts = DenseAxisArray(filter(:tec => ==("ts"), pTec).eff_st, sPeer)
    par_eff_ch_ts = DenseAxisArray(filter(:tec => ==("ts"), pTec).eff_ch, sPeer)    
    JuMP.@constraint(m, ecOpt_ts[iPeer ∈ sPeer, iY ∈ sY, iTS ∈ sTS],
        + vOpt_ts_soc[iPeer, iY, sTS_next[iTS]]
        ==
        + vOpt_ts_soc[iPeer, iY, iTS] * par_eff_st_ts[iPeer]^(dT/hours_in_day)
        + dT*(
            + vOpt_ts_ch[iPeer, iY, iTS] * sqrt(par_eff_ch_ts[iPeer])
            - vOpt_ts_dch[iPeer, iY, iTS] / sqrt(par_eff_ch_ts[iPeer])
        )
    )
    ## CHP: electrical-thermal generation and fuel consumption
    ### CHPNG (natural gas)
    JuMP.@constraint(m, icOpt_chpng,
        vOpt_chpng_gen_el .≤ repeat(vCap[:, :, Tec("chpng")], inner=(1, 1, nTS))
    )
    JuMP.@constraint(m, ecCon_chpng,
        vOpt_chpng_gen_el .== vCon_chpng_ng .* repeat(filter(:tec => ==("chpng"), pTec).eff_el, inner=(1, nY, nTS))
    )
    if bCHPThCurtail
        JuMP.@constraint(m, ecOpt_chpng,
            vOpt_chpng_gen_th ./ repeat(filter(:tec => ==("chpng"), pTec).eff_th, inner=(1, nY, nTS))
            .≤
            vOpt_chpng_gen_el ./ repeat(filter(:tec => ==("chpng"), pTec).eff_el, inner=(1, nY, nTS))
        )
    else
        JuMP.@constraint(m, ecOpt_chpng,
            vOpt_chpng_gen_th ./ repeat(filter(:tec => ==("chpng"), pTec).eff_th, inner=(1, nY, nTS))
            .==
            vOpt_chpng_gen_el ./ repeat(filter(:tec => ==("chpng"), pTec).eff_el, inner=(1, nY, nTS))
        )
    end
    ### CHPH2 (hydrogen)
        JuMP.@constraint(m, icOpt_chph2,
        vOpt_chph2_gen_el .≤ repeat(vCap[:, :, Tec("chph2")], inner=(1, 1, nTS))
    )
    JuMP.@constraint(m, ecCon_chph2,
        vOpt_chph2_gen_el .== vCon_chph2_h2 .* repeat(filter(:tec => ==("chph2"), pTec).eff_el, inner=(1, nY, nTS))
    )
    if bCHPThCurtail
        JuMP.@constraint(m, ecOpt_chph2,
            vOpt_chph2_gen_th ./ repeat(filter(:tec => ==("chph2"), pTec).eff_th, inner=(1, nY, nTS))
            .≤
            vOpt_chph2_gen_el ./ repeat(filter(:tec => ==("chph2"), pTec).eff_el, inner=(1, nY, nTS))
        )
    else
        JuMP.@constraint(m, ecOpt_chph2,
            vOpt_chph2_gen_th ./ repeat(filter(:tec => ==("chph2"), pTec).eff_th, inner=(1, nY, nTS))
            .==
            vOpt_chph2_gen_el ./ repeat(filter(:tec => ==("chph2"), pTec).eff_el, inner=(1, nY, nTS))
        )
    end
    # Demand shift ----------------------------------------------------------------
    par_dem_el_ev_peak = DenseAxisArray(dropdims(maximum(reshape(pYTS.dem_el_ev, nPeer, nY, nTS), dims=3), dims=3), sPeer, sY)
    par_dem_el_ev_mean = DenseAxisArray(dropdims(sum(reshape(pYTS.dem_el_ev, nPeer, nY, nTS) / nTS, dims=3), dims=3), sPeer, sY)
    for iPeer ∈ sPeer, iY ∈ sY
        limit_increase = par_dem_el_ev_peak[iPeer, iY] * gp(pSca, :value, [iPeer.value, "flexev_share"])
        limit_soc = par_dem_el_ev_mean[iPeer, iY] * dT * gp(pSca, :value, [iPeer.value, "flexev_share"]) * gp(pSca, :value, [iPeer.value, "flexev_duration"])
        ## Maximum increase or decrease as a share of peak of planned demand
        JuMP.set_lower_bound.(vDem_evshift_increase[iPeer, iY, :], max.(-limit_increase, -m[:dem_el_ev][iPeer, iY, :]))
        JuMP.set_upper_bound.(vDem_evshift_increase[iPeer, iY, :], +limit_increase)
        ## Maximum and minimum of shifted energy (SOC) as a share * shift duration of mean of planned demand
        JuMP.set_lower_bound.(vDem_evshift_soc[iPeer, iY, :], -limit_soc)
        JuMP.set_upper_bound.(vDem_evshift_soc[iPeer, iY, :], +limit_soc)
    end
    JuMP.@constraint(m, ecDem_evshift[iPeer ∈ sPeer, iY ∈ sY, iTS ∈ sTS],
        vDem_evshift_soc[iPeer, iY, sTS_next[iTS]]
        ==
        vDem_evshift_soc[iPeer, iY, iTS] + dT*vDem_evshift_increase[iPeer, iY, iTS]
    )
    # Physical energy balance -----------------------------------------------------
    ## Electricity balance
    ### Peer level: physical balance
    JuMP.@constraint(m, ecBalph_peer_el,
        .+ vXCph_peer_elImp .- vXCph_peer_elExp
        .+ vFS_source_el .- vFS_sink_el
        .==
        .+ m[:dem_el_hh]
        .+ m[:dem_el_ev] .+ vDem_evshift_increase
        .- vOpt_pv_gen .- vOpt_wt_gen .- vOpt_chpng_gen_el .- vOpt_chph2_gen_el
        .+ vCon_hp_el
        .+ vOpt_es_ch .- vOpt_es_dch
    )
    ### Communal level
    JuMP.@constraint(m, ecBalph_comm_el[iY ∈ sY, iTS ∈ sTS],
        + vXCph_comm_elImp[iY, iTS] - vXCph_comm_elExp[iY, iTS]
        ==
        + sum(vXCph_peer_elImp[:, iY, iTS]) - sum(vXCph_peer_elExp[:, iY, iTS])
    )
    ## Thermal balance
    ### Peer level: physical balance
    JuMP.@constraint(m, ecBalph_peer_th,
        .+ vXCph_peer_thImp .- vXCph_peer_thExp
        .+ vFS_source_th .- vFS_sink_th
        .==
        .+ m[:dem_th_hh] .- vDem_thCur
        .- vOpt_gb_gen .- vOpt_hb_gen .- vOpt_hp_gen .- vOpt_chpng_gen_th .- vOpt_chph2_gen_th
        .+ vOpt_ts_ch .- vOpt_ts_dch
    )
    ### Communal level
    JuMP.@constraint(m, ecBalph_comm_th[iY ∈ sY, iTS ∈ sTS],
        + vXCph_comm_thImp[iY, iTS] - vXCph_comm_thExp[iY, iTS]
        ==
        + sum(vXCph_peer_thImp[:, iY, iTS]) - sum(vXCph_peer_thExp[:, iY, iTS])
    )
    ## Hydrogen balance
    ### Communal level
    JuMP.@constraint(m, ecBalph_comm_h2[iY ∈ sY, iTS ∈ sTS],
        + vXCph_comm_h2Imp[iY, iTS]
        ==
        + sum(vCon_hb_h2[:, iY, iTS]) + sum(vCon_chph2_h2[:, iY, iTS])
    )
    ## Natural gas balance
    ### Communal level
    JuMP.@constraint(m, ecBalph_comm_ng[iY ∈ sY, iTS ∈ sTS],
        + vXCph_comm_ngImp[iY, iTS]
        ==
        + sum(vCon_gb_ng[:, iY, iTS]) + sum(vCon_chpng_ng[:, iY, iTS])
    )
    # Accounting energy balance ---------------------------------------------------
    if bIntTrade
        JuMP.@variable(m, 0.0 ≤ vXCac_inttrade_elImp[sPeer, sY, sTS])
        JuMP.@variable(m, 0.0 ≤ vXCac_inttrade_elExp[sPeer, sY, sTS])
        JuMP.@variable(m, 0.0 ≤ vXCac_inttrade_thImp[sPeer, sY, sTS])
        JuMP.@variable(m, 0.0 ≤ vXCac_inttrade_thExp[sPeer, sY, sTS])
        JuMP.@variable(m, 0.0 ≤ vXCac_exttrade_elImp[sPeer, sY, sTS])
        JuMP.@variable(m, 0.0 ≤ vXCac_exttrade_elExp[sPeer, sY, sTS])
        JuMP.@variable(m, 0.0 ≤ vXCac_exttrade_thImp[sPeer, sY, sTS])
        JuMP.@variable(m, 0.0 ≤ vXCac_exttrade_thExp[sPeer, sY, sTS])
        ## Electricity balance
        ### Peer level: accounting balance
        JuMP.@constraint(m, ecBalac_elImp,
            .+ vXCph_peer_elImp .== .+ vXCac_inttrade_elImp .+ vXCac_exttrade_elImp
        )
        JuMP.@constraint(m, ecBalac_elExp,
            .+ vXCph_peer_elExp .== .+ vXCac_inttrade_elExp .+ vXCac_exttrade_elExp
        )
        ### Market clearing
        JuMP.@constraint(m, ecBalac_market_el[iY ∈ sY, iTS ∈ sTS],
            + sum(vXCac_inttrade_elExp[:, iY, iTS]) - sum(vXCac_inttrade_elImp[:, iY, iTS]) == 0.0
        )
        ## Thermal balance
        ### Peer level: accounting balance
        JuMP.@constraint(m, ecBalac_thImp,
            .+ vXCph_peer_thImp .== .+ vXCac_inttrade_thImp .+ vXCac_exttrade_thImp
        )
        JuMP.@constraint(m, ecBalac_thExp,
            .+ vXCph_peer_thExp .== .+ vXCac_inttrade_thExp .+ vXCac_exttrade_thExp
        )
        ### Market clearing
        JuMP.@constraint(m, ecBalac_market_th[iY ∈ sY, iTS ∈ sTS],
            + sum(vXCac_inttrade_thExp[:, iY, iTS]) - sum(vXCac_inttrade_thImp[:, iY, iTS]) == 0.0
        )
        ## Block certain peers from participating in the internal trade
        JuMP.fix.(vXCac_inttrade_elImp[sExcludedPeer, sY, sTS], 0.0, force=true)
        JuMP.fix.(vXCac_inttrade_elExp[sExcludedPeer, sY, sTS], 0.0, force=true)
        JuMP.fix.(vXCac_inttrade_thImp[sExcludedPeer, sY, sTS], 0.0, force=true)
        JuMP.fix.(vXCac_inttrade_thExp[sExcludedPeer, sY, sTS], 0.0, force=true)
    end
    # Cost expressions ------------------------------------------------------------
    if bIntTrade
        JuMP.@expression(m, cXC_comm[iY ∈ sY], 0.0)
        JuMP.@expression(m, cXC_inttrade[iPeer ∈ sPeer, iY ∈ sY],
            # Electricity trade among peers
            + m[:inttrade_fee_el][iPeer] * dT * (
                + sum(vXCac_inttrade_elImp[iPeer, iY, :])
                + sum(vXCac_inttrade_elExp[iPeer, iY, :])
            )
            # Thermal energy trade among peers
            + m[:inttrade_fee_th][iPeer] * dT * (
                + sum(vXCac_inttrade_thImp[iPeer, iY, :])
                + sum(vXCac_inttrade_thExp[iPeer, iY, :])
            )
        )
        JuMP.@expression(m, cXC_exttrade[iPeer ∈ sPeer, iY ∈ sY],
            # Electricity trade with the utility
            + m[:impprice_el][iPeer, iY] * dT * sum(vXCac_exttrade_elImp[iPeer, iY, :])
            - m[:expprice_el][iPeer, iY] * dT * sum(vXCac_exttrade_elExp[iPeer, iY, :])
            # Thermal energy trade with the utility
            + m[:impprice_th][iPeer, iY] * dT * sum(vXCac_exttrade_thImp[iPeer, iY, :])
            - m[:expprice_th][iPeer, iY] * dT * sum(vXCac_exttrade_thExp[iPeer, iY, :])
            # Natural gas import from the utility
            + m[:impprice_ng][iPeer, iY] * dT * (
                + sum(vCon_gb_ng[iPeer, iY, :]) + sum(vCon_chpng_ng[iPeer, iY, :])
            )
            # Hydrogen import from the utility
            + m[:impprice_h2][iPeer, iY] * dT * (
                + sum(vCon_hb_h2[iPeer, iY, :]) + sum(vCon_chph2_h2[iPeer, iY, :])
            )
        )
    else
        JuMP.@expression(m, cXC_comm[iY ∈ sY],
            + minimum(m[:impprice_el][:, iY]) * dT * sum(vXCph_comm_elImp[iY, :])
            + minimum(m[:impprice_th][:, iY]) * dT * sum(vXCph_comm_thImp[iY, :])
            + minimum(m[:impprice_ng][:, iY]) * dT * sum(vXCph_comm_ngImp[iY, :])
            + minimum(m[:impprice_h2][:, iY]) * dT * sum(vXCph_comm_h2Imp[iY, :])
            - maximum(m[:expprice_el][:, iY]) * dT * sum(vXCph_comm_elExp[iY, :])
            - maximum(m[:expprice_th][:, iY]) * dT * sum(vXCph_comm_thExp[iY, :])
        )
        JuMP.@expression(m, cXC_inttrade[iPeer ∈ sPeer, iY ∈ sY], 0.0)
        JuMP.@expression(m, cXC_exttrade[iPeer ∈ sPeer, iY ∈ sY], 0.0)
    end
    par_scap_finalpaymentyear = DenseAxisArray(reshape(pTec.scap_finalpaymentyear, nPeer, nTec), sPeer, sTec)
    par_scap_payment = DenseAxisArray(reshape(pTec.scap_payment, nPeer, nTec), sPeer, sTec)
    par_ncap_loanduration = DenseAxisArray(reshape(pTec.ncap_loanduration, nPeer, nTec), sPeer, sTec)
    par_ncap_annuityrate = DenseAxisArray(reshape(pTec.ncap_annuityrate, nPeer, nTec), sPeer, sTec)
    JuMP.@expression(m, cCAPEX[iPeer ∈ sPeer, iY ∈ sY],
        + sum(
            (
                + sum(par_ncap_annuityrate[iPeer, iTec] * vInv[iPeer, iiY, iTec] for iiY ∈ sY if (iiY ≤ iY) && (iY < iiY + par_ncap_loanduration[iPeer, iTec]))
                + ((iY ≤ par_scap_finalpaymentyear[iPeer, iTec]) ? par_scap_payment[iPeer, iTec] : 0.0)
            ) for iTec ∈ sTec
        ) * (nTS/timesteps_in_year) # Reduce annual lump sum costs by the considered time horizon
    )
    par_cost_dec = DenseAxisArray(reshape(pTec.cost_dec, nPeer, nTec), sPeer, sTec)
    JuMP.@expression(m, cDec_scap[iPeer ∈ sPeer, iY ∈ sY],
        + sum(
            (
                + zDec_scap[iPeer, iY, iTec] * par_cost_dec[iPeer, iTec]
            ) for iTec ∈ sTec
        ) * (nTS/timesteps_in_year) # Reduce annual lump sum costs by the considered time horizon
    )
    par_cost_var = DenseAxisArray(reshape(pTec.cost_var, nPeer, nTec), sPeer, sTec)
    par_cost_fix = DenseAxisArray(reshape(pTec.cost_fix, nPeer, nTec), sPeer, sTec)
    JuMP.@expression(m, cOPEX[iPeer ∈ sPeer, iY ∈ sY],
        # Fixed operation costs
        + sum(
            (
                + vCap[iPeer, iY, iTec] * par_cost_fix[iPeer, iTec]
            ) for iTec ∈ sTec
        ) * (nTS/timesteps_in_year) # Reduce annual lump sum costs by the considered time horizon
        # Variable operation costs without fuel costs
        + dT * (
                + par_cost_var[iPeer, Tec("pv")] * sum(vOpt_pv_gen[iPeer, iY, :])
                + par_cost_var[iPeer, Tec("wt")] * sum(vOpt_wt_gen[iPeer, iY, :])
                + par_cost_var[iPeer, Tec("es")] * sum(vOpt_es_dch[iPeer, iY, :])
                + par_cost_var[iPeer, Tec("gb")] * sum(vOpt_gb_gen[iPeer, iY, :])
                + par_cost_var[iPeer, Tec("hb")] * sum(vOpt_hb_gen[iPeer, iY, :])
                + par_cost_var[iPeer, Tec("hp")] * sum(vOpt_hp_gen[iPeer, iY, :])
                + par_cost_var[iPeer, Tec("dh")] * sum(vXCph_peer_thImp[iPeer, iY, :])
                + par_cost_var[iPeer, Tec("ts")] * sum(vOpt_ts_dch[iPeer, iY, :])
                + par_cost_var[iPeer, Tec("chpng")] * sum(vOpt_chpng_gen_el[iPeer, iY, :])
                + par_cost_var[iPeer, Tec("chph2")] * sum(vOpt_chph2_gen_el[iPeer, iY, :])
        )
    )
    JuMP.@expression(m, cFS[iPeer ∈ sPeer, iY ∈ sY],
        + dT * bigM * sum(vFS_sink_el[iPeer, iY, :])
        + dT * bigM * sum(vFS_source_el[iPeer, iY, :])
        + dT * bigM * sum(vFS_sink_th[iPeer, iY, :])
        + dT * bigM * sum(vFS_source_th[iPeer, iY, :])
    )
    # Additional configuration ----------------------------------------------------
    ## Oneoff investment
    if bOneoff
        JuMP.@variable(m, zInv_oneoff[sPeer, sY, sTec], Bin)
        JuMP.@constraint(m, icOneoff1, vCap_new .≤ bigM .* zInv_oneoff)
        JuMP.@constraint(m, icOneoff2[iPeer ∈ sPeer, iTec ∈ sTec],
            sum(zInv_oneoff[iPeer, iY, iTec] for iY ∈ sY) ≤ 1.0
        )
    end
    ## Enforce the communal fuelswitch
    if bFuelswitch
        JuMP.@variable(m, zXC_comm_ngImp[sY], Bin)
        JuMP.@variable(m, zXC_comm_h2Imp[sY], Bin)
        JuMP.@constraint(m, icXC_comm_ngImp[iY ∈ sY],
            dT*sum(vXCph_comm_ngImp[iY, :]) ≤ bigM * nTS * zXC_comm_ngImp[iY]
        )
        JuMP.@constraint(m, icXC_comm_h2Imp[iY ∈ sY],
            dT*sum(vXCph_comm_h2Imp[iY, :]) ≤ bigM * nTS * zXC_comm_h2Imp[iY]
        )
        JuMP.@constraint(m, icXC_comm_fuelswitch[iY ∈ sY], zXC_comm_ngImp[iY] + zXC_comm_h2Imp[iY] ≤ 1.0)
    end
    ## Thermal demand curtailment
    if bConElas
        JuMP.@constraint(m, icDem_thCur_limit,
            vDem_thCur_rate .≤ repeat(filter(:parameter => ==("thCur_max"), pSca).value, inner=(1, nY))
        )
        JuMP.@constraint(m, ecDem_thCur,
            vDem_thCur .== reshape(pYTS.dem_th_hh, nPeer, nY, nTS) .* repeat(vDem_thCur_rate, inner=(1, 1, nTS))
        )
        par_dem_th_hh_total = dT .* DenseAxisArray(sum(reshape(pYTS.dem_th_hh, nPeer, nY, nTS), dims=3)[:, :, 1], sPeer, sY)
        JuMP.@expression(m, cElas[iPeer ∈ sPeer, iY ∈ sY],
            + gp(pSca, :value, [iPeer.value, "thCur_baseprice"]) * (vDem_thCur_rate[iPeer, iY] * par_dem_th_hh_total[iPeer, iY])
            # The exponent must be an integer; otherwise, JuMP will not recognize the expression as quadratic.
            + gp(pSca, :value, [iPeer.value, "thCur_incrementprice"])/2.0 * (vDem_thCur_rate[iPeer, iY] * par_dem_th_hh_total[iPeer, iY])^2
        )
    else
        JuMP.fix.(vDem_thCur, 0.0, force=true)
        JuMP.fix.(vDem_thCur_rate, 0.0, force=true)
        JuMP.@expression(m, cElas[iPeer ∈ sPeer, iY ∈ sY], 0.0)
    end
    ## No expansion
    if bNoExpand
        JuMP.delete.(m, icDec_scap); JuMP.unregister.(m, :icDec_scap)
        JuMP.fix.(vCap_new, 0.0, force=true)
        JuMP.fix.(zDec_scap, 0.0, force=true)
    end
    ## Internal trading
    if bIntTrade
        JuMP.unset_binary.(zDec_scap)
        JuMP.fix.(zDec_scap, 0.0, force=true)
    end
    # Objective -------------------------------------------------------------------
    JuMP.@objective(m, Min,
        + sum(cXC_comm)
        + sum(cXC_inttrade) + sum(cXC_exttrade)
        + sum(cCAPEX)
        + sum(cOPEX)
        + sum(cDec_scap)
        + sum(cElas)
        + sum(cFS)
    )
    return m
end
