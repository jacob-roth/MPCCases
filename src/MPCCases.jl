module MPCCases
using DelimitedFiles, StructArrays
## based on: https://github.com/StructJuMP/StructJuMP.jl/tree/master/examples/PowerGrid

include("load.jl")
export load_case
export mapLinesToBuses
export mapBusIdToIdx, mapIdxToBusId
export mapGenersToBuses
export computeAdmitances
export OPFData
export CaseData

end # module
