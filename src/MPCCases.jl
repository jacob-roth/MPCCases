module MPCCases
using DelimitedFiles
## based on: https://github.com/StructJuMP/StructJuMP.jl/tree/master/examples/PowerGrid

include("load.jl")
export loadcase
export mapLinesToBuses
export mapBusIdToIdx
export mapGenersToBuses

end # module
