module MPCCases
using DelimitedFiles, StructArrays
using SparseArrays
using Random
## based on: https://github.com/StructJuMP/StructJuMP.jl/tree/master/examples/PowerGrid

include("load.jl")
export load_case
export mapLinesToBuses
export mapBusIdToIdx, mapIdxToBusId
export mapGenersToBuses
export computeAdmitances, computeAdmittances, computeAdmittanceMatrix
export OPFData
export CaseData

include("adjust.jl")
export adj_params, adj_multi_params

include("map_id.jl")
export apply_mapping_to_file

include("util.jl")
export complete_file_path, fill_write_file_path, match_length, cp_remaining_files, get_phys
export get_Bs_psd_adjustments, update_case!

end # module
