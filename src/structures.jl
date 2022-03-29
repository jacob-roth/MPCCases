struct FileStructure
    read_file_path::String
    file_name::String
    file_ext::String
    write_file_path::String
    write_file_name::String
end

struct Adjustments
    column_names::Union{String, AbstractArray{String}}
    prod_fac::Real
    add_fac::Real
    discard_neg_vals::Bool
end

struct GaussianNoise
    mean::Union{Nothing,Real}
    sd::Union{Nothing,Real}
    rng::Union{Nothing,AbstractRNG}
end

abstract type AbstractFileIndices end

struct BusFileIndices <: AbstractFileIndices
    P::Int
    Q::Int
end

struct GenFileIndices <: AbstractFileIndices
    P::Int
    Q::Int
end

struct GenCostFileIndices <: AbstractFileIndices
    c2::Int
    c1::Int
    c0::Int
end

struct BranchFileIndices <: AbstractFileIndices
    rateA::Int
end

BusFileIndices() = BusFileIndices(3,4)
GenFileIndices() = GenFileIndices(2,3)
GenCostFileIndices() = GenCostFileIndices(5,6,7)
BranchFileIndices() = BranchFileIndices(6)
