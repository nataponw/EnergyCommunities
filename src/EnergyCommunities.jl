module EnergyCommunities

# Import dependent package ====================================================
import JuMP, JuMP.Containers.DenseAxisArray
import DataFrames
import HDF5, SQLite, DBInterface

# Export lists ================================================================
# Define set structs
export Year, Peer, Tec
# process databases and parameters
export create_blankdatabase, process_parameters
# Initialize a model
export initialize_model
# Load and unload model results
export save_results, load_results
# Process cost allocation
export iterate_over_all_coalition, allocate_cost_contribution, allocate_cost_marginalprice
# Miscellaneous functions
export convert_modelobj2df

include(joinpath(@__DIR__, "SetStructs.jl"))

# Constants ===================================================================
const bigM = 100.0
const sTec = Tec.(sort(["pv", "wt", "es", "gb", "hb", "hp", "dh", "ts", "chpng", "chph2"]))
const hours_in_year = 8760
const hours_in_day = 24.0

# Include source code =========================================================
include(joinpath(@__DIR__, "ProcessParameters.jl"))
include(joinpath(@__DIR__, "ModelInitialize.jl"))
include(joinpath(@__DIR__, "ModelIO.jl"))
include(joinpath(@__DIR__, "CostAllocation.jl"))
include(joinpath(@__DIR__, "Miscellaneous.jl"))

end
