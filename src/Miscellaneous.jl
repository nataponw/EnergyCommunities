"""
    convert_modelobj2df(obj)

Convert a model-related denseaxisarray into a dataframe with formated column names
"""
function convert_modelobj2df(obj)
    dims = list_dim_denseaxisarray(obj)
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

"""
    gp(df::DataFrames.DataFrame, selcol::Symbol, filter::Vector)

Get a parameter (gp) from a dataframe `df` based on a filter
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
    list_dim_denseaxisarray(obj)

List dimensions of a denseaxisarray object
"""
function list_dim_denseaxisarray(obj)
    dims = String[]
    for i ∈ 1:ndims(obj)
        index = first(obj.axes[i])
        if isa(index, Peer)
            push!(dims, "sPeer")
        elseif isa(index, Year)
            push!(dims, "sY")
        elseif isa(index, Int)
            push!(dims, "sTS")
        elseif isa(index, Tec)
            push!(dims, "sTec")
        end
    end
    return dims
end
