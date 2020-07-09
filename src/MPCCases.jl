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
export OPFData
export CaseData

include("admittance.jl")
export computeAdmitances, computeAdmittances, computeAdmittanceMatrix

include("adjust.jl")
export adj_params, adj_multi_params

include("compose.jl")
export compose_file

include("paste_vals.jl")
export get_paste_vals, get_shed_bus_idx

include("map_id.jl")
export apply_mapping_to_file

include("util.jl")
export complete_file_path, complete_base_files
export get_y_idx, get_write_cols_idx, fill_write_file_path, match_length, cp_remaining_files, get_path_dict
export add_gaussian_noise, undo_neg_vals
export get_phys
export get_Bs_psd_adjustments, update_case!

end # module
