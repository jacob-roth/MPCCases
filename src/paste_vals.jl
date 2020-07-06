# Generate Paste Values for adjust.jl

function get_paste_vals(arr::VecOrMat{<:Real}, file_ext::String, P::Bool, Q::Bool, x_idx::Union{Nothing, Vector{Int}, Colon}, val_to_paste::Union{Int, Float64})
    if isnothing(x_idx)
        write_cols_idx = get_write_cols_idx(file_ext)
        return arr[:, write_cols_idx]
    else
        y_idx = get_y_idx(file_ext, P, Q)
        arr[x_idx, y_idx] .= val_to_paste
        return arr[:, y_idx]
    end
end

function get_paste_vals(read_file_path::String, file_name::String, file_ext::String, P::Bool, Q::Bool, x_idx::Union{Nothing, Vector{Int}, Colon}, val_to_paste::Union{Int, Float64}; T::Type=Float64)
    read_file_path = complete_file_path(read_file_path)
    arr = readdlm(read_file_path * file_name * file_ext, T)
    return get_paste_vals(arr, file_ext, P, Q, x_idx, val_to_paste)
end

function get_paste_vals(arr::VecOrMat{<:Real}, file_ext::String, c2::Bool, c1::Bool, c0::Bool, x_idx::Union{Nothing, Vector{Int}, Colon}, val_to_paste::Union{Int, Float64})
    if isnothing(x_idx)
        write_cols_idx = get_write_cols_idx(file_ext)
        return arr[:, write_cols_idx]
    else
        y_idx = get_y_idx(file_ext, c2, c1, c0)
        arr[x_idx, y_idx] .= val_to_paste
        return arr[:, y_idx]
    end
end

function get_paste_vals(read_file_path::String, file_name::String, file_ext::String, c2::Bool, c1::Bool, c0::Bool, x_idx::Union{Nothing, Vector{Int}, Colon}, val_to_paste::Union{Int, Float64}; T::Type=Float64)
    read_file_path = complete_file_path(read_file_path)
    arr = readdlm(read_file_path * file_name * file_ext, T)
    return get_paste_vals(arr, file_ext, c2, c1, c0, x_idx, val_to_paste)
end

function get_paste_vals(arr::VecOrMat{<:Real}, file_ext::String, rateA::Bool, x_idx::Union{Nothing, Vector{Int}, Colon}, val_to_paste::Union{Int, Float64})
    if isnothing(x_idx)
        write_cols_idx = get_write_cols_idx(file_ext)
        return arr[:, write_cols_idx]
    else
        y_idx = get_y_idx(file_ext, rateA)
        arr[x_idx, y_idx] .= val_to_paste
        return arr[:, y_idx]
    end
end

function get_paste_vals(read_file_path::String, file_name::String, file_ext::String, rateA::Bool, x_idx::Union{Nothing, Vector{Int}, Colon}, val_to_paste::Union{Int, Float64}; T::Type=Float64)
    read_file_path = complete_file_path(read_file_path)
    arr = readdlm(read_file_path * file_name * file_ext, T)
    return get_paste_vals(arr, file_ext, rateA, x_idx, val_to_paste)
end

## Complementary Functions

function get_path_dict(cascades_root::String, case_name::String)
    cascades_root = complete_file_path(cascades_root)
    path_dict = Dict{Symbol, String}()
    path_dict[:casedata_path] = cascades_root * "casedata/" * case_name * "/"
    path_dict[:operatingdata_path] = cascades_root * "operatingdata/" * case_name * "/"
    path_dict[:kmcdata_path] = cascades_root * "kmcdata/"
    return path_dict
end

function get_rates(cascades_root::String, case_name::String, oppt_dir_name::String)
    path_dict = get_path_dict(cascades_root, case_name)
    rates = vec(readdlm(path_dict[:operatingdata_path] * oppt_dir_name * "/rates.csv"))
    return rates
end

function get_case_file(cascades_root::String, case_name::String, file_name::String, file_ext::String)
    path_dict = get_path_dict(cascades_root, case_name)
    case_file = readdlm(path_dict[:casedata_path] * file_name * file_ext)
    return case_file
end

function get_shed_line_idx(cascades_root::String, case_name::String, oppt_dir_name::String; 
                           rate_thresh::Union{Nothing, Int, Float64}=nothing, num_lines::Union{Nothing, Int}=nothing)
    rates = get_rates(cascades_root, case_name, oppt_dir_name)
    sorted_rate_idx = sortperm(rates, rev=true)
    if !isnothing(rate_thresh) & isnothing(num_lines)
        shed_line_idx = sorted_rate_idx[rates[sorted_rate_idx] .>= rate_thresh]
    elseif isnothing(rate_thresh) & !isnothing(num_lines)
        shed_line_idx = sorted_rate_idx[1:num_lines]
    elseif !isnothing(rate_thresh) & !isnothing(num_lines)
        shed_line_idx_tmp = sorted_rate_idx[rates[sorted_rate_idx] .>= rate_thresh]
        num_lines_tmp = min(length(shed_line_idx_tmp), num_lines)
        shed_line_idx = shed_line_idx_tmp[1:num_lines_tmp]
    else
        shed_line_idx = sorted_rate_idx
    end
    return shed_line_idx
end

function get_shed_bus_idx(cascades_root::String, case_name::String, oppt_dir_name::String, file_name::String; 
                          rate_thresh::Union{Nothing, Int, Float64}=nothing, num_lines::Union{Nothing, Int}=nothing)
    shed_line_idx = get_shed_line_idx(cascades_root, case_name, oppt_dir_name, rate_thresh=rate_thresh, num_lines=num_lines)
    branches = get_case_file(cascades_root, case_name, file_name, ".branch")
    shed_bus_idx = convert(Array{Int}, unique(branches[shed_line_idx, [1,2]]))
    return shed_bus_idx
end
