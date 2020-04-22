# Adjusting Parameters

function adj_params(read_file_path::String, file_name::String, file_ext::String,  P::Bool, Q::Bool, prod_fac::Union{Int, Float64}=1, add_fac::Union{Int, Float64}=0; start_x_idx::Int=1, end_x_idx::Int=0, T::Type=Float64, mean::Union{Nothing, Int, Float64}=nothing, sd::Union{Nothing, Int, Float64}=nothing, overwrite_file::Bool=false, write_file_path::String="", seed::Union{Nothing, Int}=nothing)
    read_file_path = complete_file_path(read_file_path)
    arr = readdlm(read_file_path * file_name * file_ext, T)
    y_idx = get_y_idx(file_ext, P, Q)
    end_x_idx = end_x_idx >= start_x_idx ? end_x_idx : size(arr, 1)
    vals = generate_vals(arr, start_x_idx, end_x_idx, y_idx, prod_fac, add_fac)
    perturbed_vals = isnothing(mean) | isnothing(sd) ? vals : add_gaussian_noise(vals, mean, sd, seed)
    adj_arr = adj_vals_in_arr(arr, start_x_idx, y_idx, perturbed_vals)
    filled_write_file_path = fill_write_file_path(write_file_path, read_file_path, overwrite_file, "adj/")
    open(filled_write_file_path * file_name * file_ext, "w") do io
        writedlm(io, arr)
    end
end

function adj_params(read_file_path::String, file_name::String, file_ext::String,  P::Bool, Q::Bool, vals::VecOrMat{<:Real}=zeros(Int,0); start_x_idx::Int=1,  T::Type=Float64, mean::Union{Nothing, Int, Float64}=nothing, sd::Union{Nothing, Int, Float64}=nothing, overwrite_file::Bool=false, write_file_path::String="", seed::Union{Nothing, Int}=nothing)
    read_file_path = complete_file_path(read_file_path)
    arr = readdlm(read_file_path * file_name * file_ext, T)
    y_idx = get_y_idx(file_ext, P, Q)
    vals = isempty(vals) ? arr[start_x_idx:end, y_idx] : reshape_vals(vals, P, Q)
    perturbed_vals = isnothing(mean) | isnothing(sd) ? vals : add_gaussian_noise(vals, mean, sd, seed)
    adj_arr = adj_vals_in_arr(arr, start_x_idx, y_idx, perturbed_vals)
    filled_write_file_path = fill_write_file_path(write_file_path, read_file_path, overwrite_file, "adj/")
    open(filled_write_file_path * file_name * file_ext, "w") do io
        writedlm(io, arr)
    end
end

function adj_params(read_file_path::String, file_name::String, file_ext::String, c2::Bool, c1::Bool, c0::Bool, prod_fac::Union{Int, Float64}=1, add_fac::Union{Int, Float64}=0; start_x_idx::Int=1, end_x_idx::Int=0, T::Type=Float64, mean::Union{Nothing, Int, Float64}=nothing, sd::Union{Nothing, Int, Float64}=nothing, overwrite_file::Bool=false, write_file_path::String="", seed::Union{Nothing, Int}=nothing)
    read_file_path = complete_file_path(read_file_path)
    arr = readdlm(read_file_path * file_name * file_ext, T)
    y_idx = get_y_idx(file_ext, c2, c1, c0)
    end_x_idx = end_x_idx >= start_x_idx ? end_x_idx : size(arr, 1)
    vals = generate_vals(arr, start_x_idx, end_x_idx, y_idx, prod_fac, add_fac)
    perturbed_vals = isnothing(mean) | isnothing(sd) ? vals : add_gaussian_noise(vals, mean, sd, seed)
    adj_arr = adj_vals_in_arr(arr, start_x_idx, y_idx, perturbed_vals)
    filled_write_file_path = fill_write_file_path(write_file_path, read_file_path, overwrite_file, "adj/")
    open(filled_write_file_path * file_name * file_ext, "w") do io
        writedlm(io, arr)
    end
end

function adj_params(read_file_path::String, file_name::String, file_ext::String, c2::Bool, c1::Bool, c0::Bool, vals::VecOrMat{<:Real}=zeros(Int,0); start_x_idx::Int=1, T::Type=Float64, mean::Union{Nothing, Int, Float64}=nothing, sd::Union{Nothing, Int, Float64}=nothing, overwrite_file::Bool=false, write_file_path::String="", seed::Union{Nothing, Int}=nothing)
    read_file_path = complete_file_path(read_file_path)
    arr = readdlm(read_file_path * file_name * file_ext, T)
    y_idx = get_y_idx(file_ext, c2, c1, c0)
    vals = isempty(vals) ? arr[start_x_idx:end, y_idx] : reshape_vals(vals, c2, c1, c0)
    perturbed_vals = isnothing(mean) | isnothing(sd) ? vals : add_gaussian_noise(vals, mean, sd, seed)
    adj_arr = adj_vals_in_arr(arr, start_x_idx, y_idx, perturbed_vals)
    filled_write_file_path = fill_write_file_path(write_file_path, read_file_path, overwrite_file, "adj/")
    open(filled_write_file_path * file_name * file_ext, "w") do io
        writedlm(io, arr)
    end
end

function adj_params(read_file_path::String, file_name::String, file_ext::Union{String, NTuple{N, String}}, P::Union{Bool, Tuple{Bool, Vararg{Bool}}}, Q::Union{Bool, Tuple{Bool, Vararg{Bool}}}, c2::Bool, c1::Bool, c0::Bool, prod_fac::Union{Int, Float64, NTuple{N, <:Union{Int, Float64}}}=1, add_fac::Union{Int, Float64, NTuple{N, <:Union{Int, Float64}}}=0; start_x_idx::Int=1, end_x_idx::Int=0, T::Type=Float64, mean::Union{Nothing, Int, Float64, NTuple{N, <:Union{Int, Float64}}}=nothing, sd::Union{Nothing, Int, Float64, NTuple{N, <:Union{Int, Float64}}}=nothing, overwrite_file::Bool=false, write_file_path::String="", seed::Union{Nothing, Int}=nothing) where {N}
    if isa(file_ext, Tuple{Vararg{String}})
        @assert file_ext == Tuple(unique(file_ext))
        PQ_file_ext = filter(x -> x in [".bus", ".gen"], [file_ext...])
        for idx in 1:length(file_ext)
            f_ext = file_ext[idx]
            p_fac = isa(prod_fac, Union{Int, Float64}) ? prod_fac : prod_fac[idx]
            a_fac = isa(add_fac, Union{Int, Float64}) ? add_fac : add_fac[idx]
            m = isnothing(mean) | isa(mean, Union{Int, Float64}) ? mean : mean[idx]
            s = isnothing(sd) | isa(sd, Union{Int, Float64}) ? sd : sd[idx]
            if f_ext in [".bus", ".gen"]
                PQ_idx = findfirst(x -> x == f_ext, PQ_file_ext)
                p, q = P[PQ_idx], Q[PQ_idx]
                adj_params(read_file_path, file_name, f_ext, p, q, p_fac, a_fac, start_x_idx=start_x_idx, end_x_idx=end_x_idx, T=T, mean=m, sd=s, overwrite_file=overwrite_file, write_file_path=write_file_path, seed=seed)
            elseif f_ext == ".gencost"
                adj_params(read_file_path, file_name, f_ext, c2, c1, c0, p_fac, a_fac, start_x_idx=start_x_idx, end_x_idx=end_x_idx, T=T, mean=m, sd=s, overwrite_file=overwrite_file, write_file_path=write_file_path, seed=seed)
            else
                throw(DomainError("file_ext is not properly defined."))
            end
        end
    else
        @assert isa(P, Bool) & isa(Q, Bool) & isa(prod_fac, Union{Int, Float64}) & isa(add_fac, Union{Int, Float64}) & isa(mean, Union{Nothing, Int, Float64}) & isa(sd, Union{Nothing, Int, Float64})
        if file_ext in [".bus", ".gen"]
            adj_params(read_file_path, file_name, file_ext, P, Q, prod_fac, add_fac, start_x_idx=start_x_idx, T=T, mean=mean, sd=sd, overwrite_file=overwrite_file, write_file_path=write_file_path, seed=seed)
        elseif file_ext == ".gencost"
            adj_params(read_file_path, file_name, file_ext, c2, c1, c0, prod_fac, add_fac, start_x_idx=start_x_idx, T=T, mean=mean, sd=sd, overwrite_file=overwrite_file, write_file_path=write_file_path, seed=seed)
        else
            throw(DomainError("file_ext is not properly defined."))
        end
    end
end

function adj_params(read_file_path::String, file_name::String, file_ext::Union{String, NTuple{N, String}}, P::Union{Bool, Tuple{Bool, Vararg{Bool}}}, Q::Union{Bool, Tuple{Bool, Vararg{Bool}}}, c2::Bool, c1::Bool, c0::Bool, vals::Union{VecOrMat{<:Real}, NTuple{N, VecOrMat{<:Real}}}; start_x_idx::Int=1, T::Type=Float64, mean::Union{Nothing, Int, Float64, NTuple{N, <:Union{Int, Float64}}}=nothing, sd::Union{Nothing, Int, Float64, NTuple{N, <:Union{Int, Float64}}}=nothing, overwrite_file::Bool=false, write_file_path::String="", seed::Union{Nothing, Int}=nothing) where {N}
    if isa(file_ext, Tuple{Vararg{String}})
        @assert file_ext == Tuple(unique(file_ext))
        PQ_file_ext = filter(x -> x in [".bus", ".gen"], [file_ext...])
        for idx in 1:length(file_ext)
            f_ext = file_ext[idx]
            m = isnothing(mean) | isa(mean, Union{Int, Float64}) ? mean : mean[idx]
            s = isnothing(sd) | isa(sd, Union{Int, Float64}) ? sd : sd[idx]
            val = isa(vals, VecOrMat{<:Real}) ? vals : vals[idx]
            if f_ext in [".bus", ".gen"]
                PQ_idx = findfirst(x -> x == f_ext, PQ_file_ext)
                p, q = P[PQ_idx], Q[PQ_idx]
                adj_params(read_file_path, file_name, f_ext, p, q, val, start_x_idx=start_x_idx, T=T, mean=m, sd=s, overwrite_file=overwrite_file, write_file_path=write_file_path, seed=seed)
            elseif f_ext == ".gencost"
                adj_params(read_file_path, file_name, f_ext, c2, c1, c0, val, start_x_idx=start_x_idx, T=T, mean=m, sd=s, overwrite_file=overwrite_file, write_file_path=write_file_path, seed=seed)
            else
                throw(DomainError("file_ext is not properly defined."))
            end
        end
    else
        @assert isa(P, Bool) & isa(Q, Bool) & isa(vals, VecOrMat{<:Real}) & isa(mean, Union{Nothing, Int, Float64}) & isa(sd, Union{Nothing, Int, Float64})
        if file_ext in [".bus", ".gen"]
            adj_params(read_file_path, file_name, file_ext, P, Q, vals, start_x_idx=start_x_idx, T=T, mean=mean, sd=sd, overwrite_file=overwrite_file, write_file_path=write_file_path, seed=seed)
        elseif file_ext == ".gencost"
            adj_params(read_file_path, file_name, file_ext, c2, c1, c0, vals, start_x_idx=start_x_idx, T=T, mean=mean, sd=sd, overwrite_file=overwrite_file, write_file_path=write_file_path, seed=seed)
        else
            throw(DomainError("file_ext is not properly defined."))
        end
    end
end

# Helper Functions for Adjusting Parameters

function complete_file_path(file_path::String)
    return !(file_path[end] ∈ Set(['/',"/"])) ? file_path * "/" : file_path
end

function fill_write_file_path(curr_write_file_path::String, read_file_path::String, overwrite_file::Bool, suffix::String)
    filled_write_file_path = complete_file_path(mkpath(
        overwrite_file                  ?   read_file_path  :
        !isempty(curr_write_file_path)  ?   curr_write_file_path :
                                            read_file_path * suffix))
    return filled_write_file_path
end

function get_y_idx(file_ext::String, P::Bool, Q::Bool)
    if file_ext == ".bus"
        return P & Q ? (3:4) : P ? (3:3) : (4:4)
    elseif file_ext == ".gen"
        return P & Q ? (2:3) : P ? (2:2) : (3:3)
    else
        throw(DomainError(file_ext, "file_ext is not properly defined for the number of boolean parameters given."))
    end
end

function get_y_idx(file_ext::String, c2::Bool, c1::Bool, c0::Bool)
    if file_ext == ".gencost"
        return c2 & c1 & c0 ? (5:7) : c2 & c1 ? (5:6) : c2 & c0 ? (5:2:7) : c1 & c0 ? (6:7) : c2 ? (5:5) : c1 ? (6:6) : (7:7)
    else
        throw(DomainError(file_ext, "file_ext is not properly defined for the number of boolean parameters given."))
    end
end

function generate_vals(arr::Array{<:Real, 2}, start_x_idx::Int, end_x_idx::Int, y_idx::OrdinalRange{<:Real}, prod_fac::Union{Int,Float64}, add_fac::Union{Int, Float64})
    subset_arr = arr[start_x_idx:end_x_idx, y_idx]
    return (prod_fac .* subset_arr) .+ add_fac
end

function reshape_vals(vals::Array{<:Real}, P::Bool, Q::Bool)
    if size(vals, 2) == P + Q
        return vals
    elseif size(vals, 2) == 1
        if P & Q
            return reshape(vals, :, 2)
        elseif P | Q
            return reshape(vals, :, 1)
        else
            throw(DomainError("vals is not properly defined for the number of boolean parameters given."))
        end
    else
        throw(DomainError("vals is not properly defined for the number of boolean parameters given."))
    end
end

function reshape_vals(vals::Array{<:Real}, c2::Bool, c1::Bool, c0::Bool)
    if size(vals, 2) == c2 + c1 + c0
        return vals
    elseif size(vals, 2) == 1
        if c2 & c1 & c0
            return reshape(vals, :, 3)
        elseif (c2 & c1) | (c2 & c0) | (c1 & c0)
            return reshape(vals, :, 2)
        elseif c2 | c1 | c0
            return reshape(vals, :, 1)
        else
            throw(DomainError("vals is not properly defined for the number of boolean parameters given."))
        end
    else
        throw(DomainError("vals is not properly defined for the number of boolean parameters given."))
    end
end

function add_gaussian_noise(vals::Array, mean::Union{Int, Float64}, sd::Union{Int, Float64}, seed::Union{Nothing, Int})
    rng = isnothing(seed) ? MersenneTwister() : MersenneTwister(seed)
    dims = size(vals)
    gaussian_noise = randn(rng, dims)
    scaled_gaussian_noise = (sd .* gaussian_noise) .+ mean
    return vals + scaled_gaussian_noise
end

function adj_vals_in_arr(arr::Array{<:Real,2}, start_x_idx::Int, y_idx::OrdinalRange{<:Real}, vals::Array)
    vals_length = size(vals, 1)
    arr[start_x_idx : start_x_idx+vals_length-1, y_idx] = vals
    return arr
end

# Adjusting Parameters and Writing Back Out to Disk

function adj_params_and_write(adj_case_name::String, adj_case_path::String, adj_case_ext::String, P::Bool, Q::Bool, prod_fac::Union{Int, Float64}, add_fac::Union{Int, Float64}; start_x_idx::Int=1, end_x_idx::Int=0, T::Type=Float64, mean::Union{Nothing, Int, Float64}=nothing, sd::Union{Nothing, Int, Float64}=nothing, overwrite_file::Bool=false, write_file_path::String="", seed::Union{Nothing, Int}=nothing)
    adj_params(adj_case_path, adj_case_name, adj_case_ext, P, Q, prod_fac, add_fac, start_x_idx=start_x_idx, end_x_idx=end_x_idx, T=T, mean=mean, sd=sd, overwrite_file=overwrite_file, write_file_path=write_file_path, seed=seed)
    cp_remaining_files(adj_case_path, write_file_path, adj_case_name)
end

function adj_params_and_write(adj_case_name::String, adj_case_path::String, adj_case_ext::String, P::Bool, Q::Bool, vals::VecOrMat{<:Real}; start_x_idx::Int=1, T::Type=Float64, mean::Union{Nothing, Int, Float64}=nothing, sd::Union{Nothing, Int, Float64}=nothing, overwrite_file::Bool=false, write_file_path::String="", seed::Union{Nothing, Int}=nothing)
    adj_params(adj_case_path, adj_case_name, adj_case_ext, P, Q, vals, start_x_idx=start_x_idx, T=T, mean=mean, sd=sd, overwrite_file=overwrite_file, write_file_path=write_file_path, seed=seed)
    cp_remaining_files(adj_case_path, write_file_path, adj_case_name)
end

function adj_params_and_write(adj_case_name::String, adj_case_path::String, adj_case_ext::String, c2::Bool, c1::Bool, c0::Bool, prod_fac::Union{Int, Float64}, add_fac::Union{Int, Float64}; start_x_idx::Int=1, end_x_idx::Int=0, T::Type=Float64, mean::Union{Nothing, Int, Float64}=nothing, sd::Union{Nothing, Int, Float64}=nothing, overwrite_file::Bool=false, write_file_path::String="", seed::Union{Nothing, Int}=nothing)
    adj_params(adj_case_path, adj_case_name, adj_case_ext, c2, c1, c0, prod_fac, add_fac, start_x_idx=start_x_idx, end_x_idx=end_x_idx, T=T, mean=mean, sd=sd, overwrite_file=overwrite_file, write_file_path=write_file_path, seed=seed)
    cp_remaining_files(adj_case_path, write_file_path, adj_case_name)
end

function adj_params_and_write(adj_case_name::String, adj_case_path::String, adj_case_ext::String, c2::Bool, c1::Bool, c0::Bool, vals::VecOrMat{<:Real}; start_x_idx::Int=1, T::Type=Float64, mean::Union{Nothing, Int, Float64}=nothing, sd::Union{Nothing, Int, Float64}=nothing, overwrite_file::Bool=false, write_file_path::String="", seed::Union{Nothing, Int}=nothing)
    adj_params(adj_case_path, adj_case_name, adj_case_ext, c2, c1, c0, vals, start_x_idx=start_x_idx, T=T, mean=mean, sd=sd, overwrite_file=overwrite_file, write_file_path=write_file_path, seed=seed)
    cp_remaining_files(adj_case_path, write_file_path, adj_case_name)
end

function adj_params_and_write(adj_case_name::String, adj_case_path::String, adj_case_ext::NTuple{N, String}, P::Union{Bool, Tuple{Bool, Vararg{Bool}}}, Q::Union{Bool, Tuple{Bool, Vararg{Bool}}}, c2::Bool, c1::Bool, c0::Bool, prod_fac::Union{Int, Float64, NTuple{N, <:Union{Int, Float64}}}, add_fac::Union{Int, Float64, NTuple{N, <:Union{Int, Float64}}}; start_x_idx::Int=1, end_x_idx::Int=0, T::Type=Float64, mean::Union{Nothing, Int, Float64, NTuple{N, <:Union{Int, Float64}}}=nothing, sd::Union{Nothing, Int, Float64, NTuple{N, <:Union{Int, Float64}}}=nothing, overwrite_file::Bool=false, write_file_path::String="", seed::Union{Nothing, Int}=nothing) where {N}
    adj_params(adj_case_path, adj_case_name, adj_case_ext, P, Q, c2, c1, c0, prod_fac, add_fac, start_x_idx=start_x_idx, end_x_idx=end_x_idx, T=T, mean=mean, sd=sd, overwrite_file=overwrite_file, write_file_path=write_file_path, seed=seed)
    cp_remaining_files(adj_case_path, write_file_path, adj_case_name)
end

function adj_params_and_write(adj_case_name::String, adj_case_path::String, adj_case_ext::NTuple{N, String}, P::Union{Bool, Tuple{Bool, Vararg{Bool}}}, Q::Union{Bool, Tuple{Bool, Vararg{Bool}}}, c2::Bool, c1::Bool, c0::Bool, vals::NTuple{N, VecOrMat{<:Real}}; start_x_idx::Int=1, T::Type=Float64, mean::Union{Nothing, Int, Float64, NTuple{N, <:Union{Int, Float64}}}=nothing, sd::Union{Nothing, Int, Float64, NTuple{N, <:Union{Int, Float64}}}=nothing, overwrite_file::Bool=false, write_file_path::String="", seed::Union{Nothing, Int}=nothing) where {N}
    adj_params(adj_case_path, adj_case_name, adj_case_ext, P, Q, c2, c1, c0, vals, start_x_idx=start_x_idx, T=T, mean=mean, sd=sd, overwrite_file=overwrite_file, write_file_path=write_file_path, seed=seed)
    cp_remaining_files(adj_case_path, write_file_path, adj_case_name)
end

# Helper Function for Adjusting Parameters and Writing Back to Disk

function cp_remaining_files(src_path::String, dst_path::String, file_name::String)
    file_exts = [".bus", ".gen", ".gencost", ".branch", ".phys"]
    src_file_path = complete_file_path(src_path) * file_name
    dst_file_path = complete_file_path(dst_path) * file_name
    for ext in file_exts
        if !isfile(dst_file_path * ext)
            cp(src_file_path * ext, dst_file_path * ext, force=false)
        end
    end
end
