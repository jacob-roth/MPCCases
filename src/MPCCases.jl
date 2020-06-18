module MPCCases
using DelimitedFiles, StructArrays
using SparseArrays
using Random
## based on: https://github.com/StructJuMP/StructJuMP.jl/tree/master/examples/PowerGrid

include("load.jl")
export load_case, complete_base_files
export mapLinesToBuses
export mapBusIdToIdx, mapIdxToBusId
export mapGenersToBuses
export OPFData
export CaseData

include("admittance.jl")
export computeAdmitances, computeAdmittances, computeAdmittanceMatrix

include("adjust.jl")
export adj_params, adj_multi_params

include("compose.jl")
export compose_file

include("map_id.jl")
export apply_mapping_to_file

include("util.jl")
export complete_file_path, get_y_idx, get_write_cols_idx, fill_write_file_path, match_length, cp_remaining_files, get_phys
export get_Bs_psd_adjustments, update_case!

end # module
