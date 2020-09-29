# Generate Paste Values for adjust.jl

# Shed Load based on OPF solution

# 1: P and Q
function get_paste_vals(arr::VecOrMat{<:Real}, file_ext::String, P::Bool, Q::Bool, x_idx::Union{Nothing, Vector{Int}, Colon}, 
                        val_to_paste::Union{Int, Float64})
    if isnothing(x_idx)
        write_cols_idx = get_write_cols_idx(file_ext)
        return arr[:, write_cols_idx]
    else
        y_idx = get_y_idx(file_ext, P, Q)
        arr[x_idx, y_idx] .= val_to_paste
        return arr[:, y_idx]
    end
end

# 2: P and Q
function get_paste_vals(read_file_path::String, file_name::String, file_ext::String, P::Bool, Q::Bool, x_idx::Union{Nothing, Vector{Int}, Colon}, 
                        val_to_paste::Union{Int, Float64}; 
                        T::Type=Float64)
    read_file_path = complete_file_path(read_file_path)
    arr = readdlm(read_file_path * file_name * file_ext, T)
    return get_paste_vals(arr, file_ext, P, Q, x_idx, val_to_paste)
end

# 3: c2, c1, c0
function get_paste_vals(arr::VecOrMat{<:Real}, file_ext::String, c2::Bool, c1::Bool, c0::Bool, x_idx::Union{Nothing, Vector{Int}, Colon}, 
                        val_to_paste::Union{Int, Float64})
    if isnothing(x_idx)
        write_cols_idx = get_write_cols_idx(file_ext)
        return arr[:, write_cols_idx]
    else
        y_idx = get_y_idx(file_ext, c2, c1, c0)
        arr[x_idx, y_idx] .= val_to_paste
        return arr[:, y_idx]
    end
end

# 4: c2, c1, c0
function get_paste_vals(read_file_path::String, file_name::String, file_ext::String, c2::Bool, c1::Bool, c0::Bool, x_idx::Union{Nothing, Vector{Int}, Colon}, 
                        val_to_paste::Union{Int, Float64}; 
                        T::Type=Float64)
    read_file_path = complete_file_path(read_file_path)
    arr = readdlm(read_file_path * file_name * file_ext, T)
    return get_paste_vals(arr, file_ext, c2, c1, c0, x_idx, val_to_paste)
end

# 5: rateA
function get_paste_vals(arr::VecOrMat{<:Real}, file_ext::String, rateA::Bool, x_idx::Union{Nothing, Vector{Int}, Colon}, 
                        val_to_paste::Union{Int, Float64})
    if isnothing(x_idx)
        write_cols_idx = get_write_cols_idx(file_ext)
        return arr[:, write_cols_idx]
    else
        y_idx = get_y_idx(file_ext, rateA)
        arr[x_idx, y_idx] .= val_to_paste
        return arr[:, y_idx]
    end
end

# 6: rateA
function get_paste_vals(read_file_path::String, file_name::String, file_ext::String, rateA::Bool, x_idx::Union{Nothing, Vector{Int}, Colon}, 
                        val_to_paste::Union{Int, Float64}; 
                        T::Type=Float64)
    read_file_path = complete_file_path(read_file_path)
    arr = readdlm(read_file_path * file_name * file_ext, T)
    return get_paste_vals(arr, file_ext, rateA, x_idx, val_to_paste)
end

# Shed Load for Initial Contingencies

function get_paste_vals(read_file_path::String, file_name::String, file_ext::String, P::Bool, Q::Bool, x_idx::Union{Nothing, Vector{Int}, Colon}, 
                        prod_fac::Union{Int, Float64}, add_fac::Union{Int, Float64};
                        T::Type=Float64)
    read_file_path = complete_file_path(read_file_path)
    arr = readdlm(read_file_path * file_name * file_ext, T)
    y_idx = get_y_idx(file_ext, P, Q)
    if isnothing(x_idx)
        write_cols_idx = get_write_cols_idx(file_ext)
        return arr[:, write_cols_idx]
    else
        arr[x_idx, y_idx] = prod_fac * arr[x_idx, y_idx] .+ add_fac
        return arr[:, y_idx]
    end
end

## Complementary Functions

function get_shed_line_idx(cascades_root::String, case_name::String, oppt_dir_name::String; 
                           rate_thresh::Union{Nothing, Int, Float64}=nothing, num_lines::Union{Nothing, Int}=nothing)
    rates = get_xbar_arr(cascades_root, case_name, oppt_dir_name, "rates")
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
    branches = get_case_arr(cascades_root, case_name, file_name, ".branch")
    shed_bus_idx = convert(Array{Int}, unique(branches[shed_line_idx, [1,2]]))
    return shed_bus_idx
end
