function adj_params(file_structure::FileStructure, adjustments::Adjustments, gaussian_noise::GaussianNoise)
    arr = readdlm(complete_file_path(file_structure.read_file_path) * file_structure.file_name * file_structure.file_ext)
    y_idx = get_y_idx(file_structure.file_ext, adjustments.column_names)
    vals = generate_vals(arr, y_idx, adjustments.prod_fac, adjustments.add_fac)
    perturbed_vals = isnothing(gaussian_noise.mean) | isnothing(gaussian_noise.sd) ? vals : 
                     add_gaussian_noise(vals, gaussian_noise.mean, gaussian_noise.sd, gaussian_noise.rng)
    arr[:, y_idx] = perturbed_vals
    abs_arr = adjustments.discard_neg_vals ? undo_neg_vals(arr, y_idx, vals) : arr
    filled_write_file_path = fill_write_file_path(file_structure.write_file_path, complete_file_path(file_structure.read_file_path), "adj/")
    filled_write_file_name = fill_write_file_name(file_structure.file_name, file_structure.write_file_name)
    write_cols_idx = get_write_cols_idx(file_structure.file_ext)
    open(filled_write_file_path * filled_write_file_name * file_structure.file_ext, "w") do io
        writedlm(io, abs_arr[:, write_cols_idx])
    end
end

# Helper Functions for Adjusting Parameters
function get_y_idx(file_ext::String, column_name::String)
    if file_ext == ".bus"
        bus_file_indices = BusFileIndices()
        if column_name == "P"
            return bus_file_indices.P
        elseif column_name == "Q"
            return bus_file_indices.Q
        else
            throw(DomainError(file_ext, "column_name not found in file_ext."))
        end
    elseif file_ext == ".gen"
        gen_file_indices = GenFileIndices()
        if column_name == "P"
            return gen_file_indices.P
        elseif column_name == "Q"
            return gen_file_indices.Q
        else
            throw(DomainError(file_ext, "column_name not found in file_ext."))
        end
    elseif file_ext == ".gencost"
        gencost_file_indices = GenCostFileIndices()
        if column_name == "c2"
            return gencost_file_indices.c2
        elseif column_name == "c1"
            return gencost_file_indices.c1
        elseif column_name == "c0"
            return gencost_file_indices.c0
        else 
            throw(DomainError(file_ext, "column_name not found in file_ext."))
        end
    elseif file_ext == ".branch"
        branch_file_indices = BranchFileIndices()
        if column_name == "rateA"
            return branch_file_indices.rateA
        else
            throw(DomainError(file_ext, "column_name not found in file_ext."))
        end
    else
        throw(DomainError(file_ext, "file_ext not defined."))
    end
end

function get_y_idx(file_ext::String, column_names::AbstractArray{String})
    y_idx = [get_y_idx(file_ext, column_name) for column_name in column_names]
    return y_idx
end

function generate_vals(arr::AbstractArray, y_idx::Union{Real, AbstractArray}, prod_fac::Real, add_fac::Real)
    subset_arr = arr[:, y_idx]
    return (prod_fac .* subset_arr) .+ add_fac
end

function add_gaussian_noise(vals::AbstractArray, mean::Real, sd::Real, rng::Union{Nothing, AbstractRNG})
    dims = size(vals)
    gaussian_noise = randn(rng, dims)
    scaled_gaussian_noise = (sd .* gaussian_noise) .+ mean
    return vals + scaled_gaussian_noise
end

# For Pd, c2, c1, c0, if adjustment in adj_arr is negative, use the values from vals instead. Also used in Pg from xbar
function undo_neg_vals(arr::AbstractArray, y_idx::Union{Real, AbstractArray}, vals::AbstractArray)
    subset_arr = arr[:, y_idx]
    discard_neg_vals = subset_arr .* (subset_arr .>= 0) + vals .* (subset_arr .< 0)
    arr[:, y_idx] = discard_neg_vals
    return arr
end

function get_write_cols_idx(file_ext::String)
    if file_ext ∈ (".bus", ".gen")
        column_names = ["P", "Q"]
    elseif file_ext == ".gencost"
        column_names = ["c2", "c1", "c0"]
    elseif file_ext == ".branch"
        column_names = "rateA"
    else
        throw(DomainError(file_ext, "file_ext is not properly defined."))
    end
    return get_y_idx(file_ext, column_names)
end
