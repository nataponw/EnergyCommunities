"""
    create_blankdatabase(sY_int::Vector{Int}, nTS::Int; dbname="blank.db")

Create a black database based on `db_structure_peer.sql` and populate with mandatory keys
"""
function create_blankdatabase(sY_int::Vector{Int}, nTS::Int; dbname="blank.db")
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
    process_parameters(dbPaths::Vector{String}, sY::Vector{Year}, sDay::Vector{Int}, dT::Real)

Extract and filter parameters from all peers listed in `dbPaths`

# Arguments
- `dbPaths` : Paths to DB of individual peers
- `sY` : Selected years
- `sDay` : Selected days
- `dT` : Time interval used to map `sTS`

See also : [`process_parametersingleAgent`](@ref)
"""
function process_parameters(dbPaths::Vector{String}, sY::Vector{Year}, sDay::Vector{Int}, dT::Real)
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
        par_peer = process_parameters_singlepeer(dbpath, sY, sTS)
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

"""
    extractDBTable(db, tableName::String)

Extract a table `tableName` from a DB connection `db` as a DataFrame
"""
function extractDBTable(db, tableName::String)
    df = DataFrames.DataFrame(DBInterface.execute(db, "SELECT * FROM ($tableName)"))
    return df
end

"""
    process_parameters_singlepeer(dbpath::String, sY::Vector{Year}, sTS::Vector{Int})

Extract and filter parameter tables from `dbpath` for individual peers
"""
function process_parameters_singlepeer(dbpath::String, sY::Vector{Year}, sTS::Vector{Int})
    db = SQLite.DB(dbpath)
    agentid = extractDBTable(db, "id")
    pSca = extractDBTable(db, "pSca")
    pY = extractDBTable(db, "pY")
    filter!(:year => ∈(sY), pY)
    pTec = extractDBTable(db, "pTec")
    filter!(:tec => ∈(sTec), pTec)
    pYTec = extractDBTable(db, "pYTec")
    filter!([:year, :tec] => (y, t) -> (y ∈ sY) & (t ∈ sTec), pYTec)
    pYTS = extractDBTable(db, "pYTS")
    filter!([:year, :timestep] => (y, ts) -> (y ∈ sY) & (ts ∈ sTS), pYTS)
    return (; agentid, pSca, pY, pTec, pYTec, pYTS)
end
