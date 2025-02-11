"""
    save_results(m::JuMP.Model, filename::String; bSaveConstraintDual::Bool=false)

Save all results and essential sets into a HDF5 file

# Keyword Arguments
- `bSaveConstraintDual` : save dual variables

See also : [`load_results`](@ref)
"""
function save_results(m::JuMP.Model, filename::String; bSaveConstraintDual::Bool=false)
    function _write_value_dim(objname, value, dims; fid=fid)
        gid = HDF5.create_group(fid, string(objname))
        HDF5.write_dataset(gid, "value", value)
        HDF5.write_dataset(gid, "dims", dims)
        HDF5.close(gid)
        return nothing
    end
    # Create or overwrite an existing file
    fid = HDF5.h5open(filename, "w")
    # Categorize registered objects
    vArrVar = Symbol[]; vArrCon = Symbol[]; vArrExp = Symbol[]; vArrPar = Symbol[]
    vVar = Symbol[]; vCon = Symbol[]; vExp = Symbol[]; vPar = Symbol[]; vSet = Symbol[]
    for obj ∈ keys(JuMP.object_dictionary(m))
        objType = typeof(m[obj])
        if objType <: JuMP.Containers.DenseAxisArray
            elementType = typeof(first(m[obj]))
            if elementType <: JuMP.VariableRef
                push!(vArrVar, obj)
            elseif elementType <: JuMP.ConstraintRef
                push!(vArrCon, obj)
            elseif elementType <: JuMP.AffExpr
                push!(vArrExp, obj)
            elseif elementType <: Number
                push!(vArrPar, obj)
            else
                @warn "Unable to identify the type of $(string(obj))"
            end
        elseif objType <: JuMP.VariableRef
            push!(vVar, obj)
        elseif objType <: JuMP.ConstraintRef
            push!(vCon, obj)
        elseif objType <: JuMP.AffExpr
            push!(vExp, obj)
        elseif objType <: Number
            push!(vPar, obj)
        elseif objType <: Vector
            push!(vSet, obj)
        else
            @warn "Unable to identify the type of $(string(obj))"
        end
    end
    # Write sets
    HDF5.write_dataset(fid, "_index_sPeer", [iPeer.value for iPeer ∈ m[:sPeer]])
    HDF5.write_dataset(fid, "_index_sY", [iY.value for iY ∈ m[:sY]])
    HDF5.write_dataset(fid, "_index_sTS", m[:sTS])
    HDF5.write_dataset(fid, "_index_sTec", [iTec.value for iTec ∈ m[:sTec]])
    # write data
    [_write_value_dim(obj, JuMP.value.(m[obj]).data, list_dim_denseaxisarray(m[obj])) for obj ∈ vArrVar]
    [_write_value_dim(obj, JuMP.value.(m[obj]).data, list_dim_denseaxisarray(m[obj])) for obj ∈ vArrExp]
    [_write_value_dim(obj, m[obj].data, list_dim_denseaxisarray(m[obj])) for obj ∈ vArrPar]
    [_write_value_dim(obj, JuMP.value(m[obj]), ["scalar"]) for obj ∈ vVar]
    [_write_value_dim(obj, JuMP.value(m[obj]), ["scalar"]) for obj ∈ vExp]
    [_write_value_dim(obj, m[obj], ["scalar"]) for obj ∈ vPar]
    if bSaveConstraintDual
        [_write_value_dim(obj, JuMP.dual.(m[obj]).data, list_dim_denseaxisarray(m[obj])) for obj ∈ vArrCon]
        [_write_value_dim(obj, JuMP.dual(m[obj]), ["scalar"]) for obj ∈ vCon]
    end
    # Write objective values
    _write_value_dim("objvalue", JuMP.objective_value(m), ["scalar"])
    # Close the connection
    HDF5.close(fid)
    return nothing
end

"""
    load_results(filename::String; objlist::Vector{String})

Load all saved results from an HDF5 into a dictionary

# Keyword Arguments
- `objlist` : List of objects to be loaded. By default, all objects are loaded.

See also : [`save_results`](@ref)
"""
function load_results(filename::String; objlist::Vector{String}=String[])
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
