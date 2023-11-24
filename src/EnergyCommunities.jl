module EnergyCommunities

# Import dependent package ====================================================
import JuMP, JuMP.Containers.DenseAxisArray
import DataFrames
import HDF5, SQLite, DBInterface

# Export lists ================================================================
# Set structs
export Year, Peer, Tec
# Functions related to the extraction and processing of parameters
export createblankdatabase, processParameters
# Functions related to the initialization of an optimization problem
export initializeModel
# JuMP object conversion functions
export convert_modelobj2df
# Export and import functions for solved problems
export saveResults, loadResults
# Analysis
export iterate_over_all_coalition, calculate_contribution_chipin

# Custom structs ==============================================================
include(joinpath(@__DIR__, "SetStructs.jl"))

# Constants ===================================================================
const bigM = 100.0
const sTec = Tec.(sort(["pv", "wt", "es", "gb", "hb", "hp", "dh", "ts", "chpng", "chph2"]))
const hours_in_year = 8760
const hours_in_day = 24.0
const emissionfactor = Dict{Tuple{String, Year}, Float64}()

# Parameter extraction ========================================================
"""
    createblankdatabase(sY_int::Vector{Int}, nTS::Int; dbname="blank.db")

Create a black database based on `db_structure_peer.sql` and populate with mandatory keys
"""
function createblankdatabase(sY_int::Vector{Int}, nTS::Int; dbname="blank.db")
    sTec_str = map(x -> x.value, sTec); sY_int = sort(sY_int)
    db = SQLite.DB(dbname)
    f = open(joinpath(@__DIR__, "auxiliary_src/db_structure.sql"), "r")
    allcommand = read(f, String)
    close(f)
    allcommand = replace(allcommand, '\r' => "")
    allcommand = replace(allcommand, '\n' => "")
    allcommand_break = string.(split(allcommand, ";"))
    for sqlcommand ∈ allcommand_break
        isempty(sqlcommand) && continue
        DBInterface.execute(db, sqlcommand)
    end
    # Populate pY with years
    [DBInterface.execute(db, "INSERT INTO pY (year) VALUES ($(iY))") for iY ∈ sY_int]
    # Populate pTec with mandatory technologies
    [DBInterface.execute(db, "INSERT INTO pTec (tec) VALUES ('$(iTec)')") for iTec ∈ sTec_str]
    # Populate pYTec with years and technologies
    [DBInterface.execute(db, "INSERT INTO pYTec (year, tec) VALUES ('$(iY)', '$(iTec)')") for iY ∈ sY_int, iTec ∈ sTec_str]
    # Populate PYTS with years and timesteps
    insert_template = "INSERT INTO pYTS (year, timestep) VALUES "
    [insert_template *= "(!, $(i))," for i ∈ 1:nTS]
    insert_template = insert_template[1:(end-1)]
    for iY ∈ sY_int
        command = replace(insert_template, '!' => iY)
        DBInterface.execute(db, command)
    end
end

"""
    processParameters(dbPaths::Vector{String}, sY::Vector{Year}, sDay::Vector{Int}, dT::Real)

Extract and filter parameters from all peers listed in `dbPaths`

# Arguments
- `dbPaths` : Paths to DB of individual peers
- `sY` : Selected years
- `sDay` : Selected days
- `dT` : Time interval used to map `sTS`

See also : [`processParameterSingleAgent`](@ref)
"""
function processParameters(dbPaths::Vector{String}, sY::Vector{Year}, sDay::Vector{Int}, dT::Real)
    function add_peername_and_merge(df, df_add, name)
        df_add[!, :agentid] .= name
        DataFrames.select!(df_add, :agentid, :)
        df = vcat(df, df_add)
        return df
    end
    # Sets, counts, constants
    timesteps_in_day = Int(hours_in_day/dT)
    sTS = timesteps_in_day * repeat(sDay .- 1, inner=timesteps_in_day) .+ repeat(1:timesteps_in_day, outer=length(sDay))
    # Create blank DataFrames
    agentid = DataFrames.DataFrame(); pSca = DataFrames.DataFrame(); pY = DataFrames.DataFrame(); pTec = DataFrames.DataFrame(); pYTec = DataFrames.DataFrame(); pYTS = DataFrames.DataFrame()
    # Sequentially extract and merge
    for dbpath ∈ dbPaths
        par_peer = processParameterSinglePeer(dbpath, sY, sTS)
        peername = par_peer.agentid.agentid
        agentid = vcat(agentid, par_peer.agentid)
        pSca = add_peername_and_merge(pSca, par_peer.pSca, peername)
        pY = add_peername_and_merge(pY, par_peer.pY, peername)
        pTec = add_peername_and_merge(pTec, par_peer.pTec, peername)
        pYTec = add_peername_and_merge(pYTec, par_peer.pYTec, peername)
        pYTS = add_peername_and_merge(pYTS, par_peer.pYTS, peername)
    end
    # Peer set
    sPeer = Peer.(agentid.agentid)
    return (; sTS, sTec, sPeer, agentid, pSca, pY, pTec, pYTec, pYTS)
end

# Model initialization ========================================================
"""
    initializeModel(sPeer, sY, sTS, pSca, pY, pTec, pYTec, pYTS, dT; kwargs)

Build an optimization problem

# Arguments
- `sPeer`, `sY`, `sTS`, `sTec` : sets of peers, years, timesteps, and technologies, respectively.
- `pSca`, `pY`, `pTec`, `pYTec`, `pYTS` : Extracted parameter dataframes from [`processParameters`](@ref)
- `dT` : Time interval

# Keyword Arguments
- `bOneoff`
- `bFuelswitch`
- `bConElas`
- `bNoExpand`
- `bIntTrade`
- `sExcludedPeer`
"""
function initializeModel(
    sPeer, sY, sTS, sTec, pSca, pY, pTec, pYTec, pYTS, dT;
    bOneoff::Bool=false, bFuelswitch::Bool=false, bConElas::Bool=false, bNoExpand::Bool=false,
    bIntTrade::Bool=false, sExcludedPeer=Peer[]
    )
    # Sorting indexes and parameters
    sort!(sPeer); sort!(sY); sort!(sTS); sort!(sTec)
    sort!(pSca, [:agentid])
    sort!(pY, [:year, :agentid])
    sort!(pTec, [:tec, :agentid])
    sort!(pYTec, [:tec, :year, :agentid])
    sort!(pYTS, [:timestep, :year, :agentid])
    # Building the model
    m = JuMP.Model()
    # Essential mappings and parameters
    sTS_next = Dict(sTS .=> circshift(sTS, -1))
    timesteps_in_year = Int(hours_in_year/dT)
    nPeer = length(sPeer); nY = length(sY); nTS = length(sTS); nTec = length(sTec)
    # Precalculation --------------------------------------------------------------
    # Calculate annuity rates for new investment
    DataFrames.transform!(pTec, [:ncap_loaninterest, :ncap_loanduration] => (
        (r, n) -> r ./ (1 .- (1 .+ r).^(-n))
    ) => :ncap_annuityrate)
    # Essential constants and sets permanently mapped to the model ----------------
    m[:dT] = dT
    m[:sPeer] = sPeer; m[:sY] = sY; m[:sTS] = sTS; m[:sTec] = sTec
    # To-do: attach all parameters to the model
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
    JuMP.@constraint(m, ecOpt_chpng,
        vOpt_chpng_gen_el ./ repeat(filter(:tec => ==("chpng"), pTec).eff_el, inner=(1, nY, nTS))
        .==
        vOpt_chpng_gen_th ./ repeat(filter(:tec => ==("chpng"), pTec).eff_th, inner=(1, nY, nTS))
    )
    JuMP.@constraint(m, ecCon_chpng,
        vOpt_chpng_gen_el .== vCon_chpng_ng .* repeat(filter(:tec => ==("chpng"), pTec).eff_el, inner=(1, nY, nTS))
    )
    ### CHPH2 (hydrogen)
        JuMP.@constraint(m, icOpt_chph2,
        vOpt_chph2_gen_el .≤ repeat(vCap[:, :, Tec("chph2")], inner=(1, 1, nTS))
    )
    JuMP.@constraint(m, ecOpt_chph2,
        vOpt_chph2_gen_el ./ repeat(filter(:tec => ==("chph2"), pTec).eff_el, inner=(1, nY, nTS))
        .==
        vOpt_chph2_gen_th ./ repeat(filter(:tec => ==("chph2"), pTec).eff_th, inner=(1, nY, nTS))
    )
    JuMP.@constraint(m, ecCon_chph2,
        vOpt_chph2_gen_el .== vCon_chph2_h2 .* repeat(filter(:tec => ==("chph2"), pTec).eff_el, inner=(1, nY, nTS))
    )
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
            + gp(pSca, :value, [iPeer.value, "thCur_incrementprice"])/2.0 * (vDem_thCur_rate[iPeer, iY] * par_dem_th_hh_total[iPeer, iY])^2.0
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
        + sum(cDec_scap)
        + sum(cElas)
        + sum(cFS)
    )
    return m
end

# JuMP object conversion ======================================================
"""
    convert_modelobj2df(obj)

Convert a model-related denseaxisarray into a dataframe with formated column names
"""
function convert_modelobj2df(obj)
    dims = listdimdenseaxisarray(obj)
    df = convert_denseaxisarray2dataframe(obj)
    DataFrames.rename!(df, vcat(dims, "value"))
    if "sPeer" ∈ DataFrames.names(df)
        df.sPeer = map(x -> x.value, df.sPeer)
        DataFrames.rename!(df, :sPeer => :peer)
    end
    if "sY" ∈ DataFrames.names(df)
        df.sY = map(x -> x.value, df.sY)
        DataFrames.rename!(df, :sY => :year)
    end
    if "sTS" ∈ DataFrames.names(df)
        DataFrames.rename!(df, :sTS => :timestep)
    end
    if "sTec" ∈ DataFrames.names(df)
        df.sTec = map(x -> x.value, df.sTec)
        DataFrames.rename!(df, :sTec => :tec)
    end
    return df
end

# Model import and export =====================================================
"""
    saveResults(m::JuMP.Model, filename::String; bSaveConstraintDual::Bool=false)

Save all results and essential sets into a HDF5 file

# Keyword Arguments
- `bSaveConstraintDual` : save dual variables

See also : [`loadResults`](@ref)
"""
function saveResults(m::JuMP.Model, filename::String; bSaveConstraintDual::Bool=false)
    # Create or overwrite an existing file
    fid = HDF5.h5open(filename, "w")
    # Write sets
    allSets = [:sPeer, :sY, :sTS, :sTec]
    HDF5.write_dataset(fid, "_index_sPeer", [iPeer.value for iPeer ∈ m[:sPeer]])
    HDF5.write_dataset(fid, "_index_sY", [iY.value for iY ∈ m[:sY]])
    HDF5.write_dataset(fid, "_index_sTS", m[:sTS])
    HDF5.write_dataset(fid, "_index_sTec", [iTec.value for iTec ∈ m[:sTec]])
    # List objects
    allObjects = setdiff(keys(JuMP.object_dictionary(m)), allSets)
    allConstraintRef = [x for x ∈ allObjects if ((string(x)[1:2]=="ec") | (string(x)[1:2]=="ic"))]
    !bSaveConstraintDual && setdiff!(allObjects, allConstraintRef)
    # Process objects
    for obj ∈ allObjects
        gid = HDF5.create_group(fid, string(obj))
        if obj in allConstraintRef
            HDF5.write_dataset(gid, "value", JuMP.dual.(m[obj]).data)
        elseif typeof(m[obj]) <: Number
            HDF5.write_dataset(gid, "value", m[obj])
        else
            HDF5.write_dataset(gid, "value", JuMP.value.(m[obj]).data)
        end
        listDims = listdimdenseaxisarray(m[obj])
        if isempty(listDims)
            HDF5.write_dataset(gid, "dims", ["scalar"])
        else
            HDF5.write_dataset(gid, "dims", listDims)
        end
        HDF5.close(gid)
    end
    # Write essential scalar values
    exportScalar = Dict{String, Float64}("objvalue" => JuMP.objective_value(m))
    for k ∈ keys(exportScalar)
        gid = HDF5.create_group(fid, k)
        HDF5.write_dataset(gid, "value", exportScalar[k])
        HDF5.write_dataset(gid, "dims", ["scalar"]); HDF5.close(gid)
    end
    # Close the connection
    HDF5.close(fid)
    return nothing
end

"""
    loadResults(filename::String; objlist::Vector{String})

Load all saved results from an HDF5 into a dictionary

# Keyword Arguments
- `objlist` : List of objects to be loaded. By default, all objects are loaded.

See also : [`saveResults`](@ref)
"""
function loadResults(filename::String; objlist::Vector{String}=String[])
    fid = HDF5.h5open(filename, "r")
    isempty(objlist) && (objlist = [x for x ∈ HDF5.keys(fid) if !contains(x, "_index_")])
    results = Dict{Symbol, Any}()
    # recreate the sets
    sPeer = Peer.(HDF5.read(fid, "_index_sPeer")); results[:sPeer] = sPeer;
    sY = Year.(HDF5.read(fid, "_index_sY")); results[:sY] = sY;
    sTS = HDF5.read(fid, "_index_sTS"); results[:sTS] = sTS;
    sTec = Symbol.(HDF5.read(fid, "_index_sTec")); results[:sTec] = sTec;
    for objname ∈ objlist
        # get the dimension
        dims = HDF5.read(fid[objname], "dims")
        if dims == ["scalar"]
            # process scalar
            results[Symbol(objname)] = HDF5.read(fid[objname], "value")
        else
            # process array
            results[Symbol(objname)] = reconstruct_denseaxisarray(
                HDF5.read(fid[objname], "value"),
                dims, sPeer, sY, sTS, sTec
            )
        end
    end
    HDF5.close(fid)
    return results
end

# Analyse the cost allocation =================================================
include(joinpath(@__DIR__, "CostAllocation.jl"))

# Internal functions ==========================================================

"""
    extractDBTable(db, tableName::String)

Extract a table `tableName` from a DB connection `db` as a DataFrame
"""
function extractDBTable(db, tableName::String)
    df = DataFrames.DataFrame(DBInterface.execute(db, "SELECT * FROM ($tableName)"))
    return df
end

"""
    processParameterSinglePeer(dbpath::String, sY::Vector{Year}, sTS::Vector{Int})

Extract and filter parameter tables from `dbpath` for individual peers
"""
function processParameterSinglePeer(dbpath::String, sY::Vector{Year}, sTS::Vector{Int})
    db = SQLite.DB(dbpath)
    # Agent id
    agentid = extractDBTable(db, "id")
    # Scalar parameter [...]
    pSca = extractDBTable(db, "pSca")
    # Parameter [sY]
    pY = extractDBTable(db, "pY")
    filter!(:year => ∈(sY), pY)
    # Parameter [sTec]
    pTec = extractDBTable(db, "pTec")
    filter!(:tec => ∈(sTec), pTec)
    # Parameter [sY, sTec]
    pYTec = extractDBTable(db, "pYTec")
    filter!([:year, :tec] => (y, t) -> (y ∈ sY) & (t ∈ sTec), pYTec)
    # Parameter [sY-sTS]
    pYTS = extractDBTable(db, "pYTS")
    filter!([:year, :timestep] => (y, ts) -> (y ∈ sY) & (ts ∈ sTS), pYTS)
    return (; agentid, pSca, pY, pTec, pYTec, pYTS)
end

"""
Select a scalar parameter value from a DataFrame
"""
function gp(df::DataFrames.DataFrame, selcol::Symbol, filter::Vector)
    index = collect(1:size(df)[1])
    for i in 1:length(filter)
        if !ismissing(filter[i])
            intersect!(index, findall(df[:, i] .== filter[i]))
        end
    end
    if isempty(index)
        return missing
    else
        return df[index, selcol][1]
    end
end

"""
    convert_denseaxisarray2dataframe(obj::JuMP.Containers.DenseAxisArray)

Convert any denseaxisarray into a dataframe
"""
function convert_denseaxisarray2dataframe(obj::JuMP.Containers.DenseAxisArray)
    ndims = length(obj.axes)
    df = DataFrames.DataFrame(Symbol("dim$(1)") => obj.axes[1])
    if ndims > 1
        for i ∈ 2:ndims
            df = DataFrames.crossjoin(df, DataFrames.DataFrame(Symbol("dim$(i)") => obj.axes[i]))
        end
    end
    df.value = permutedims(JuMP.value.(obj).data, reverse(1:ndims))[:]
    return df
end

"""
    listdimdenseaxisarray(obj)

List dimensions of a denseaxisarray object
"""
function listdimdenseaxisarray(obj)
    dims = String[]
    for i ∈ 1:ndims(obj)
        index = first(obj.axes[i])
        isa(index, Peer) && push!(dims, "sPeer")
        isa(index, Year) && push!(dims, "sY")
        isa(index, Int) && push!(dims, "sTS")
        isa(index, Tec) && push!(dims, "sTec")
    end
    return dims
end

"""
    reconstruct_denseaxisarray(values, dims, sPeer, sY, sTS, sTec)

Reconstruct a denseaxisarray from a array and axes
"""
function reconstruct_denseaxisarray(values, dims, sPeer, sY, sTS, sTec)
    if dims == ["sPeer", "sY", "sTec"]
        obj = DenseAxisArray(values, sPeer, sY, sTec)
    elseif dims == ["sPeer", "sY", "sTS"]
        obj = DenseAxisArray(values, sPeer, sY, sTS)
    elseif dims == ["sPeer", "sY"]
        obj = DenseAxisArray(values, sPeer, sY)
    elseif dims == ["sPeer", "sTec"]
        obj = DenseAxisArray(values, sPeer, sTec)
    elseif dims == ["sY", "sTS"]
        obj = DenseAxisArray(values, sY, sTS)
    elseif dims == ["sPeer"]
        obj = DenseAxisArray(values, sPeer)
    elseif dims == ["sY"]
        obj = DenseAxisArray(values, sY)
    else
        @error "Extraction mode is not defined."
        obj = DenseAxisArray([0.0], ["dummydimension"])
    end
    return obj
end

end
