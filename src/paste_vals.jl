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

