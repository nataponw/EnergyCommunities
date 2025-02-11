include("../src/EnergyCommunities.jl")
using .EnergyCommunities
using HiGHS, JuMP, DataFrames
###############################################################################

# Prepare the parameters ======================================================
# Initiate the solver
selSolver = optimizer_with_attributes(HiGHS.Optimizer)
# Parameters related to the community setup
sY = EnergyCommunities.Year.([2030])
sDay = Vector(1:1:365)
dT = 1.0
dbPaths = ["example/database/db_pg$(i).db" for i ∈ 1:8]
# Extract parameters from the databases
(
    sTS, sTec, sPeer,
    agentid,
    pSca, pY, pTec, pYTec, pYTS
) = process_parameters(dbPaths, sY, sDay, dT)

# Optimize the status quo scenario (SQ) =======================================
function create_model_SQ(;sPeer=sPeer, sY=sY, sTS=sTS, sTec=sTec, pSca=pSca, pY=pY, pTec=pTec, pYTec=pYTec, pYTS=pYTS, dT=dT, selSolver=selSolver)
    pSca = deepcopy(pSca)
    # Disable the EV flexible charging
    pSca[pSca.parameter .∈ Ref(["flexev_share", "flexev_duration"]), :value] .= 0.0
    m = initialize_model(
        sPeer, sY, sTS, sTec, pSca, pY, pTec, pYTec, pYTS, dT,
        bOneoff=false, bFuelswitch=false, bConElas=false, bNoExpand=true,
        bCHPThCurtail=false, bIntTrade=true, sExcludedPeer=sPeer, solverbackend=selSolver
    )
    return m
end
model_SQ = create_model_SQ()
optimize!(model_SQ)

# Optimize the electricity trading scenario (EL) ==============================
function create_model_EL(;sPeer=sPeer, sY=sY, sTS=sTS, sTec=sTec, pSca=pSca, pY=pY, pTec=pTec, pYTec=pYTec, pYTS=pYTS, dT=dT, selSolver=selSolver)
    pSca = deepcopy(pSca)
    # Disable the EV flexible charging
    pSca[pSca.parameter .∈ Ref(["flexev_share", "flexev_duration"]), :value] .= 0.0
    pTec = deepcopy(pTec)
    # Disable the district heating
    pTec[pTec.tec .== "dh", [:pot, :scap_cap]] .= 0.
    m = initialize_model(
        sPeer, sY, sTS, sTec, pSca, pY, pTec, pYTec, pYTS, dT,
        bOneoff=false, bFuelswitch=false, bConElas=false, bNoExpand=true,
        bCHPThCurtail=false, bIntTrade=true, sExcludedPeer=Peer[], solverbackend=selSolver
    )
    return m
end
model_EL = create_model_EL()
optimize!(model_EL)

# Analyze the internal payment for the cooperative scheme =====================
dc_idx_objvalue, dc_idx_costs = iterate_over_all_coalition(model_EL)
allocatedCost_cooperative = allocate_cost_contribution(dc_idx_objvalue, dc_idx_costs, sPeer)
@show allocatedCost_cooperative

# Analyze the internal payment for the competitive scheme =====================
allocatedCost_competitive = allocate_cost_marginalprice(model_EL)
@show allocatedCost_competitive

# Unload the solved model instance into HDF4 file =============================
save_results(model_SQ, "results_SQ.h5")
save_results(model_EL, "results_EL.h5")
