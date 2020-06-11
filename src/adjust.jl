# Adjusting Parameters

# 1: Single file_ext, multi column, uniform scaling update for .bus and .gen (P,Q)
function adj_params(read_file_path::String, file_name::String, file_ext::String, P::Bool, Q::Bool, prod_fac::Real, add_fac::Real; 
                    start_x_idx::Int=1, end_x_idx::Int=0, T::Type=Float64, mean::Union{Nothing, Real}=nothing, sd::Union{Nothing, Real}=nothing, 
                    overwrite_file::Bool=false, write_file_path::String="", write_file_name::String="", only_write_changed_cols::Bool=false, 
                    seed::Union{Nothing, Int}=nothing)
    read_file_path = complete_file_path(read_file_path)
    arr = readdlm(read_file_path * file_name * file_ext, T)
    y_idx = get_y_idx(file_ext, P, Q)
    end_x_idx = end_x_idx >= start_x_idx ? end_x_idx : size(arr, 1)
    vals = generate_vals(arr, start_x_idx, end_x_idx, y_idx, prod_fac, add_fac)
    perturbed_vals = isnothing(mean) | isnothing(sd) ? vals : add_gaussian_noise(vals, mean, sd, seed)
    adj_arr = adj_vals_in_arr(arr, start_x_idx, y_idx, perturbed_vals)
    filled_write_file_path = fill_write_file_path(write_file_path, read_file_path, overwrite_file, "adj/")
    filled_write_file_name = fill_write_file_name(file_name, write_file_name, overwrite_file)
    write_cols_idx = get_write_cols_idx(file_ext)
    open(filled_write_file_path * filled_write_file_name * file_ext, "w") do io
        if only_write_changed_cols
            writedlm(io, adj_arr[:, write_cols_idx])
        else
            writedlm(io, adj_arr)
        end
    end
end

# 2: Single file_ext, multi column, block value update for .bus and .gen (P,Q)
function adj_params(read_file_path::String, file_name::String, file_ext::String,  P::Bool, Q::Bool, vals::VecOrMat{<:Real}; 
                    start_x_idx::Int=1,  T::Type=Float64, mean::Union{Nothing, Real}=nothing, sd::Union{Nothing, Real}=nothing, 
                    overwrite_file::Bool=false, write_file_path::String="", write_file_name::String="", only_write_changed_cols::Bool=false, 
                    seed::Union{Nothing, Int}=nothing)
    read_file_path = complete_file_path(read_file_path)
    arr = readdlm(read_file_path * file_name * file_ext, T)
    y_idx = get_y_idx(file_ext, P, Q)
    vals = isempty(vals) ? arr[start_x_idx:end, y_idx] : reshape_vals(vals, P, Q)
    perturbed_vals = isnothing(mean) | isnothing(sd) ? vals : add_gaussian_noise(vals, mean, sd, seed)
    adj_arr = adj_vals_in_arr(arr, start_x_idx, y_idx, perturbed_vals)
    filled_write_file_path = fill_write_file_path(write_file_path, read_file_path, overwrite_file, "adj/")
    filled_write_file_name = fill_write_file_name(file_name, write_file_name, overwrite_file)
    write_cols_idx = get_write_cols_idx(file_ext)
    open(filled_write_file_path * filled_write_file_name * file_ext, "w") do io
        if only_write_changed_cols
            writedlm(io, adj_arr[:, write_cols_idx])
        else
            writedlm(io, adj_arr)
        end
    end
end

# 3: Single file_ext, multi column, uniform scaling update for .gencost (c2,c1,c0)
function adj_params(read_file_path::String, file_name::String, file_ext::String, c2::Bool, c1::Bool, c0::Bool, prod_fac::Real, add_fac::Real; 
                    start_x_idx::Int=1, end_x_idx::Int=0, T::Type=Float64, mean::Union{Nothing, Real}=nothing, sd::Union{Nothing, Real}=nothing, 
                    overwrite_file::Bool=false, write_file_path::String="", write_file_name::String="", only_write_changed_cols::Bool=false, 
                    seed::Union{Nothing, Int}=nothing)
    read_file_path = complete_file_path(read_file_path)
    arr = readdlm(read_file_path * file_name * file_ext, T)
    y_idx = get_y_idx(file_ext, c2, c1, c0)
    end_x_idx = end_x_idx >= start_x_idx ? end_x_idx : size(arr, 1)
    vals = generate_vals(arr, start_x_idx, end_x_idx, y_idx, prod_fac, add_fac)
    perturbed_vals = isnothing(mean) | isnothing(sd) ? vals : add_gaussian_noise(vals, mean, sd, seed)
    adj_arr = adj_vals_in_arr(arr, start_x_idx, y_idx, perturbed_vals)
    filled_write_file_path = fill_write_file_path(write_file_path, read_file_path, overwrite_file, "adj/")
    filled_write_file_name = fill_write_file_name(file_name, write_file_name, overwrite_file)
    write_cols_idx = get_write_cols_idx(file_ext)
    open(filled_write_file_path * filled_write_file_name * file_ext, "w") do io
        if only_write_changed_cols
            writedlm(io, adj_arr[:, write_cols_idx])
        else
            writedlm(io, adj_arr)
        end
    end
end

# 4: Single file_ext, multi column, block value update for .gencost (c2,c1,c0)
function adj_params(read_file_path::String, file_name::String, file_ext::String, c2::Bool, c1::Bool, c0::Bool, vals::VecOrMat{<:Real}; 
                    start_x_idx::Int=1, T::Type=Float64, mean::Union{Nothing, Real}=nothing, sd::Union{Nothing, Real}=nothing, 
                    overwrite_file::Bool=false, write_file_path::String="", write_file_name::String="", only_write_changed_cols::Bool=false, 
                    seed::Union{Nothing, Int}=nothing)
    read_file_path = complete_file_path(read_file_path)
    arr = readdlm(read_file_path * file_name * file_ext, T)
    y_idx = get_y_idx(file_ext, c2, c1, c0)
    vals = isempty(vals) ? arr[start_x_idx:end, y_idx] : reshape_vals(vals, c2, c1, c0)
    perturbed_vals = isnothing(mean) | isnothing(sd) ? vals : add_gaussian_noise(vals, mean, sd, seed)
    adj_arr = adj_vals_in_arr(arr, start_x_idx, y_idx, perturbed_vals)
    filled_write_file_path = fill_write_file_path(write_file_path, read_file_path, overwrite_file, "adj/")
    filled_write_file_name = fill_write_file_name(file_name, write_file_name, overwrite_file)
    write_cols_idx = get_write_cols_idx(file_ext)
    open(filled_write_file_path * filled_write_file_name * file_ext, "w") do io
        if only_write_changed_cols
            writedlm(io, adj_arr[:, write_cols_idx])
        else
            writedlm(io, adj_arr)
        end
    end
end

# 5: Single file_ext, multi column, uniform scaling update for .branch (rateA)
function adj_params(read_file_path::String, file_name::String, file_ext::String, rateA::Bool, prod_fac::Real, add_fac::Real; 
                    start_x_idx::Int=1, end_x_idx::Int=0, T::Type=Float64, mean::Union{Nothing, Real}=nothing, sd::Union{Nothing, Real}=nothing, 
                    overwrite_file::Bool=false, write_file_path::String="", write_file_name::String="", only_write_changed_cols::Bool=false,
                    seed::Union{Nothing, Int}=nothing)
    read_file_path = complete_file_path(read_file_path)
    arr = readdlm(read_file_path * file_name * file_ext, T)
    y_idx = get_y_idx(file_ext, rateA)
    end_x_idx = end_x_idx >= start_x_idx ? end_x_idx : size(arr, 1)
    vals = generate_vals(arr, start_x_idx, end_x_idx, y_idx, prod_fac, add_fac)
    perturbed_vals = isnothing(mean) | isnothing(sd) ? vals : add_gaussian_noise(vals, mean, sd, seed)
    adj_arr = adj_vals_in_arr(arr, start_x_idx, y_idx, perturbed_vals)
    filled_write_file_path = fill_write_file_path(write_file_path, read_file_path, overwrite_file, "adj/")
    filled_write_file_name = fill_write_file_name(file_name, write_file_name, overwrite_file)
    write_cols_idx = get_write_cols_idx(file_ext)
    open(filled_write_file_path * filled_write_file_name * file_ext, "w") do io
        if only_write_changed_cols
            writedlm(io, adj_arr[:, write_cols_idx])
        else
            writedlm(io, adj_arr)
        end
    end
end

# 6: Single file_ext, multi column, block value update for .branch (rateA)
function adj_params(read_file_path::String, file_name::String, file_ext::String, rateA::Bool, vals::VecOrMat{<:Real};
                    start_x_idx::Int=1,  T::Type=Float64, mean::Union{Nothing, Real}=nothing, sd::Union{Nothing, Real}=nothing, 
                    overwrite_file::Bool=false, write_file_path::String="", write_file_name::String="", only_write_changed_cols::Bool=false, 
                    seed::Union{Nothing, Int}=nothing)
    read_file_path = complete_file_path(read_file_path)
    arr = readdlm(read_file_path * file_name * file_ext, T)
    y_idx = get_y_idx(file_ext, rateA)
    vals = isempty(vals) ? arr[start_x_idx:end, y_idx] : reshape_vals(vals, rateA)
    perturbed_vals = isnothing(mean) | isnothing(sd) ? vals : add_gaussian_noise(vals, mean, sd, seed)
    adj_arr = adj_vals_in_arr(arr, start_x_idx, y_idx, perturbed_vals)
    filled_write_file_path = fill_write_file_path(write_file_path, read_file_path, overwrite_file, "adj/")
    filled_write_file_name = fill_write_file_name(file_name, write_file_name, overwrite_file)
    write_cols_idx = get_write_cols_idx(file_ext)
    open(filled_write_file_path * filled_write_file_name * file_ext, "w") do io
        if only_write_changed_cols
            writedlm(io, adj_arr[:, write_cols_idx])
        else
            writedlm(io, adj_arr)
        end
    end
end

# 7: Multi file_ext, multi column, varied scaling update
function adj_params(read_file_path::String, file_name::String, file_ext::Union{String, Tuple{String, Vararg{String}}}, 
                    P::Union{Bool, Tuple{Bool, Vararg{Bool}}}, Q::Union{Bool, Tuple{Bool, Vararg{Bool}}}, c2::Bool, c1::Bool, c0::Bool, rateA::Bool, 
                    prod_fac::Union{Real, Tuple{Real, Vararg{Real}}}, add_fac::Union{Real, Tuple{Real, Vararg{Real}}}; 
                    start_x_idx::Int=1, end_x_idx::Int=0, T::Type=Float64, mean::Union{Nothing, Real, Tuple{Real, Vararg{Real}}}=nothing, sd::Union{Nothing, Real, Tuple{Real, Vararg{Real}}}=nothing, 
                    overwrite_file::Bool=false, write_file_path::String="", write_file_name::String="", only_write_changed_cols::Bool=false, 
                    seed::Union{Nothing, Int}=nothing)
    if isa(file_ext, Tuple{Vararg{String}})
        @assert file_ext == Tuple(unique(file_ext))
        PQ_file_ext = filter(x -> x in [".bus", ".gen"], [file_ext...])
        for idx in 1:length(file_ext)
            f_ext = file_ext[idx]
            m = isnothing(mean) | isa(mean..., Real) ? Real(mean...) : mean[idx]
            s = isnothing(sd) | isa(sd..., Real) ? Real(sd...) : sd[idx]
            p_fac = isa(prod_fac..., Real) ? Real(prod_fac...) : prod_fac[idx]
            a_fac = isa(add_fac..., Real) ? Real(add_fac...) : add_fac[idx]
            if f_ext in [".bus", ".gen"]
                PQ_idx = findfirst(x -> x == f_ext, PQ_file_ext)
                p, q = P[PQ_idx], Q[PQ_idx]
                adj_params(read_file_path, file_name, f_ext, p, q, p_fac, a_fac, 
                           start_x_idx=start_x_idx, end_x_idx=end_x_idx, T=T, mean=m, sd=s, 
                           overwrite_file=overwrite_file, write_file_path=write_file_path, write_file_name=write_file_name, 
                           only_write_changed_cols=only_write_changed_cols, seed=seed)
            elseif f_ext == ".gencost"
                adj_params(read_file_path, file_name, f_ext, c2, c1, c0, p_fac, a_fac, 
                           start_x_idx=start_x_idx, end_x_idx=end_x_idx, T=T, mean=m, sd=s, 
                           overwrite_file=overwrite_file, write_file_path=write_file_path, write_file_name=write_file_name, 
                           only_write_changed_cols=only_write_changed_cols, seed=seed)
            elseif f_ext == ".branch"
                adj_params(read_file_path, file_name, f_ext, rateA, p_fac, a_fac, 
                           start_x_idx=start_x_idx, end_x_idx=end_x_idx, T=T, mean=m, sd=s, 
                           overwrite_file=overwrite_file, write_file_path=write_file_path, write_file_name=write_file_name, 
                           only_write_changed_cols=only_write_changed_cols, seed=seed)
            else
                throw(DomainError("file_ext is not properly defined."))
            end
        end
    else
        @assert isa(P..., Bool) & isa(Q..., Bool) & isa(prod_fac..., Real) & isa(add_fac..., Real) & 
                isa(mean..., Union{Nothing, Real}) & isa(sd..., Union{Nothing, Real})
        p, q, p_fac, a_fac = Bool(P...), Bool(Q...), Real(prod_fac...), Real(add_fac...)
        m = isa(mean..., Real) ? Real(mean...) : nothing
        s = isa(sd..., Real) ? Real(sd...) : nothing
        if file_ext in [".bus", ".gen"]
            adj_params(read_file_path, file_name, file_ext, p, q, p_fac, a_fac, 
                       start_x_idx=start_x_idx, T=T, mean=m, sd=s, 
                       overwrite_file=overwrite_file, write_file_path=write_file_path, write_file_name=write_file_name, 
                       only_write_changed_cols=only_write_changed_cols, seed=seed)
        elseif file_ext == ".gencost"
            adj_params(read_file_path, file_name, file_ext, c2, c1, c0, p_fac, a_fac, 
                       start_x_idx=start_x_idx, T=T, mean=m, sd=s, 
                       overwrite_file=overwrite_file, write_file_path=write_file_path, write_file_name=write_file_name, 
                       only_write_changed_cols=only_write_changed_cols, seed=seed)
        elseif file_ext == ".branch"
            adj_params(read_file_path, file_name, file_ext, rateA, p_fac, a_fac, 
                       start_x_idx=start_x_idx, T=T, mean=m, sd=s, 
                       overwrite_file=overwrite_file, write_file_path=write_file_path, write_file_name=write_file_name, 
                       only_write_changed_cols=only_write_changed_cols, seed=seed)
        else
            throw(DomainError("file_ext is not properly defined."))
        end
    end
end

# 8: Multi file_ext, multi column, varied block value update
function adj_params(read_file_path::String, file_name::String, file_ext::Union{String, Tuple{String, Vararg{String}}}, 
                    P::Union{Bool, Tuple{Bool, Vararg{Bool}}}, Q::Union{Bool, Tuple{Bool, Vararg{Bool}}}, c2::Bool, c1::Bool, c0::Bool, rateA::Bool,
                    vals::Union{VecOrMat{<:Real}, Tuple{VecOrMat{<:Real}, Vararg{VecOrMat{<:Real}}}}; 
                    start_x_idx::Int=1, T::Type=Float64, mean::Union{Nothing, Real, Tuple{Real, Vararg{Real}}}=nothing, sd::Union{Nothing, Real, Tuple{Real, Vararg{Real}}}=nothing,
                    overwrite_file::Bool=false, write_file_path::String="", write_file_name::String="", only_write_changed_cols::Bool=false, 
                    seed::Union{Nothing, Int}=nothing)
    if isa(file_ext, Tuple{Vararg{String}})
        @assert file_ext == Tuple(unique(file_ext))
        PQ_file_ext = filter(x -> x in [".bus", ".gen"], [file_ext...])
        for idx in 1:length(file_ext)
            f_ext = file_ext[idx]
            m = isnothing(mean) | isa(mean..., Real) ? Real(mean...) : mean[idx]
            s = isnothing(sd) | isa(sd..., Real) ? Real(sd...) : sd[idx]
            val = isa(vals..., VecOrMat{Real}) ? Array(vals...) : vals[idx]
            if f_ext in [".bus", ".gen"]
                PQ_idx = findfirst(x -> x == f_ext, PQ_file_ext)
                p, q = P[PQ_idx], Q[PQ_idx]
                adj_params(read_file_path, file_name, f_ext, p, q, val, 
                           start_x_idx=start_x_idx, T=T, mean=m, sd=s, 
                           overwrite_file=overwrite_file, write_file_path=write_file_path, write_file_name=write_file_name, 
                           only_write_changed_cols=only_write_changed_cols, seed=seed)
            elseif f_ext == ".gencost"
                adj_params(read_file_path, file_name, f_ext, c2, c1, c0, val, 
                           start_x_idx=start_x_idx, T=T, mean=m, sd=s, 
                           overwrite_file=overwrite_file, write_file_path=write_file_path, write_file_name=write_file_name, 
                           only_write_changed_cols=only_write_changed_cols, seed=seed)
            elseif f_ext == ".branch"
                adj_params(read_file_path, file_name, f_ext, rateA, val, 
                           start_x_idx=start_x_idx, T=T, mean=m, sd=s, 
                           overwrite_file=overwrite_file, write_file_path=write_file_path, write_file_name=write_file_name, 
                           only_write_changed_cols=only_write_changed_cols, seed=seed)
            else
                throw(DomainError("file_ext is not properly defined."))
            end
        end
    else
        @assert isa(P..., Bool) & isa(Q..., Bool) & isa(vals..., VecOrMat{Real}) & 
                isa(mean..., Union{Nothing, Real}) & isa(sd..., Union{Nothing, Real})
        p, q, val = Bool(P...), Bool(Q...), Array(vals...)
        m = isa(mean..., Real) ? Real(mean...) : nothing
        s = isa(sd..., Real) ? Real(sd...) : nothing
        if file_ext in [".bus", ".gen"]
            adj_params(read_file_path, file_name, file_ext, p, q, val, 
                       start_x_idx=start_x_idx, T=T, mean=m, sd=s, 
                       overwrite_file=overwrite_file, write_file_path=write_file_path, write_file_name=write_file_name, 
                       only_write_changed_cols=only_write_changed_cols, seed=seed)
        elseif file_ext == ".gencost"
            adj_params(read_file_path, file_name, file_ext, c2, c1, c0, val, 
                       start_x_idx=start_x_idx, T=T, mean=m, sd=s, 
                       overwrite_file=overwrite_file, write_file_path=write_file_path, write_file_name=write_file_name, 
                       only_write_changed_cols=only_write_changed_cols, seed=seed)
        elseif file_ext == ".branch"
            adj_params(read_file_path, file_name, file_ext, rateA, val, 
                       start_x_idx=start_x_idx, T=T, mean=m, sd=s, 
                       overwrite_file=overwrite_file, write_file_path=write_file_path, write_file_name=write_file_name, 
                       only_write_changed_cols=only_write_changed_cols, seed=seed)
        else
            throw(DomainError("file_ext is not properly defined."))
        end
    end
end

# Helper Functions for Adjusting Parameters

function generate_vals(arr::VecOrMat{<:Real}, start_x_idx::Int, end_x_idx::Int, y_idx::OrdinalRange{<:Real}, prod_fac::Real, add_fac::Real)
    subset_arr = arr[start_x_idx:end_x_idx, y_idx]
    return (prod_fac .* subset_arr) .+ add_fac
end

function reshape_vals(vals::VecOrMat{<:Real}, P::Bool, Q::Bool)
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

function reshape_vals(vals::VecOrMat{<:Real}, c2::Bool, c1::Bool, c0::Bool)
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

function reshape_vals(vals::VecOrMat{<:Real}, rateA::Bool)
    if size(vals, 2) == 1 & rateA
        return vals
    elseif size(vals, 2) != 1 & rateA
        return reshape(vals, :, 1)   
    else
        throw(DomainError("vals is not properly defined for the number of boolean parameters given."))
    end
end

function add_gaussian_noise(vals::VecOrMat{<:Real}, mean::Real, sd::Real, seed::Union{Nothing, Int})
    rng = isnothing(seed) ? MersenneTwister() : MersenneTwister(seed)
    dims = size(vals)
    gaussian_noise = randn(rng, dims)
    scaled_gaussian_noise = (sd .* gaussian_noise) .+ mean
    return vals + scaled_gaussian_noise
end

function adj_vals_in_arr(arr::VecOrMat{<:Real}, start_x_idx::Int, y_idx::OrdinalRange{<:Real}, vals::VecOrMat{<:Real})
    vals_length = size(vals, 1)
    arr[start_x_idx : start_x_idx+vals_length-1, y_idx] = vals
    return arr
end

# Adjusting Parameters for Multiple Cases

# 1
function adj_multi_params(adj_case_path::Union{String, Tuple{String, Vararg{String}}}, adj_case_name::Union{String, Tuple{String, Vararg{String}}}, adj_case_ext::Union{String, Tuple{String, Vararg{String}}}, 
                          P::Bool, Q::Bool, prod_fac::Union{Real, Tuple{Real, Vararg{Real}}}, add_fac::Union{Real, Tuple{Real, Vararg{Real}}}; 
                          start_x_idx::Int=1, end_x_idx::Int=0, T::Type=Float64, mean::Union{Nothing, Real, Tuple{Real, Vararg{Real}}}=nothing, sd::Union{Nothing, Real, Tuple{Real, Vararg{Real}}}=nothing, 
                          overwrite_file::Bool=false, write_file_path::Union{String, Tuple{String, Vararg{String}}}="", write_file_name::Union{String, Tuple{String, Vararg{String}}}="", only_write_changed_cols::Bool=false, 
                          seed::Union{Nothing, Int}=nothing)
    N = lcm(length(adj_case_path), length(adj_case_name), length(adj_case_ext), 
            length(prod_fac), length(add_fac), length(mean), length(sd), 
            length(write_file_path), length(write_file_name))
    adj_case_path = isa(adj_case_path, Union{String, Tuple{String}}) ? match_length(String(adj_case_path), N) : adj_case_path
    adj_case_name = isa(adj_case_name, Union{String, Tuple{String}}) ? match_length(String(adj_case_name), N) : adj_case_name
    adj_case_ext = isa(adj_case_ext, Union{String, Tuple{String}}) ? match_length(String(adj_case_ext), N) : adj_case_ext
    write_file_path = (!overwrite_file & isa(write_file_path, Union{String, Tuple{String}})) ? match_length(String(write_file_path), N) : write_file_path
    write_file_name = (!overwrite_file & isa(write_file_name, Union{String, Tuple{String}})) ? match_length(String(write_file_name), N) : write_file_name

    p_fac = isa(prod_fac, Union{Real, Tuple{Real}}) ? match_length(prod_fac..., N) : prod_fac
    a_fac = isa(add_fac, Union{Real, Tuple{Real}}) ? match_length(add_fac..., N) : add_fac
    m = isa(mean, Union{Real, Tuple{Real}}) ? match_length(mean..., N) : mean
    s = isa(sd, Union{Real, Tuple{Real}}) ? match_length(sd..., N) : sd

    rng = isnothing(seed) ? MersenneTwister() : MersenneTwister(seed)
    sub_seeds = isnothing(seed) ? match_length(seed, N) : Tuple(rand(rng, Int, N))

    @assert length(adj_case_path) == length(adj_case_name) == length(adj_case_ext) == 
            length(write_file_path) == length(write_file_name) ==
            length(p_fac) == length(a_fac) == length(m) == length(s) == length(sub_seeds)

    for idx in 1:N
        adj_params(adj_case_path[idx], adj_case_name[idx], adj_case_ext[idx], P, Q, p_fac[idx], a_fac[idx], 
        start_x_idx=start_x_idx, end_x_idx=end_x_idx, T=T, mean=m[idx], sd=s[idx], 
        overwrite_file=overwrite_file, write_file_path=write_file_path[idx], write_file_name=write_file_name[idx], 
        only_write_changed_cols=only_write_changed_cols, seed=sub_seeds[idx])
    end
end

# 2
function adj_multi_params(adj_case_path::Union{String, Tuple{String, Vararg{String}}}, adj_case_name::Union{String, Tuple{String, Vararg{String}}}, adj_case_ext::Union{String, Tuple{String, Vararg{String}}}, 
                          P::Bool, Q::Bool, vals::Union{VecOrMat{<:Real}, NTuple{N, VecOrMat{<:Real}}}; 
                          start_x_idx::Int=1, end_x_idx::Int=0, T::Type=Float64, mean::Union{Nothing, Real, Tuple{Real, Vararg{Real}}}=nothing, sd::Union{Nothing, Real, Tuple{Real, Vararg{Real}}}=nothing, 
                          overwrite_file::Bool=false, write_file_path::Union{String, Tuple{String, Vararg{String}}}="", write_file_name::Union{String, Tuple{String, Vararg{String}}}="", only_write_changed_cols::Bool=false, 
                          seed::Union{Nothing, Int}=nothing) where {N}
    N = lcm(length(adj_case_path), length(adj_case_name), length(adj_case_ext), 
            length(vals), length(mean), length(sd), 
            length(write_file_path), length(write_file_name))
    adj_case_path = isa(adj_case_path, Union{String, Tuple{String}}) ? match_length(String(adj_case_path), N) : adj_case_path
    adj_case_name = isa(adj_case_name, Union{String, Tuple{String}}) ? match_length(String(adj_case_name), N) : adj_case_name
    adj_case_ext = isa(adj_case_ext, Union{String, Tuple{String}}) ? match_length(String(adj_case_ext), N) : adj_case_ext
    write_file_path = (!overwrite_file & isa(write_file_path, Union{String, Tuple{String}})) ? match_length(String(write_file_path), N) : write_file_path
    write_file_name = (!overwrite_file & isa(write_file_name, Union{String, Tuple{String}})) ? match_length(String(write_file_name), N) : write_file_name

    val = isa(vals, Union{VecOrMat{<:Real}, Tuple{VecOrMat{<:Real}}}) ? match_length(vals..., N) : vals
    m = isa(mean, Union{Real, Tuple{Real}}) ? match_length(mean..., N) : mean
    s = isa(sd, Union{Real, Tuple{Real}}) ? match_length(sd..., N) : sd

    rng = isnothing(seed) ? MersenneTwister() : MersenneTwister(seed)
    sub_seeds = isnothing(seed) ? match_length(seed, N) : Tuple(rand(rng, Int, N))

    @assert length(adj_case_path) == length(adj_case_name) == length(adj_case_ext) == 
            length(write_file_path) == length(write_file_name) ==
            length(val) == length(m) == length(s) == length(sub_seeds)

    for idx in 1:N
        adj_params(adj_case_path[idx], adj_case_name[idx], adj_case_ext[idx], P, Q, val[idx], 
        start_x_idx=start_x_idx, end_x_idx=end_x_idx, T=T, mean=m[idx], sd=s[idx], 
        overwrite_file=overwrite_file, write_file_path=write_file_path[idx], write_file_name=write_file_name[idx], 
        only_write_changed_cols=only_write_changed_cols, seed=sub_seeds[idx])
    end
end

# 3
function adj_multi_params(adj_case_path::Union{String, Tuple{String, Vararg{String}}}, adj_case_name::Union{String, Tuple{String, Vararg{String}}}, adj_case_ext::Union{String, Tuple{String, Vararg{String}}}, 
                          c2::Bool, c1::Bool, c0::Bool, prod_fac::Union{Real, Tuple{Real, Vararg{Real}}}, add_fac::Union{Real, Tuple{Real, Vararg{Real}}}; 
                          start_x_idx::Int=1, end_x_idx::Int=0, T::Type=Float64, mean::Union{Nothing, Real, Tuple{Real, Vararg{Real}}}=nothing, sd::Union{Nothing, Real, Tuple{Real, Vararg{Real}}}=nothing, 
                          overwrite_file::Bool=false, write_file_path::Union{String, Tuple{String, Vararg{String}}}="", write_file_name::Union{String, Tuple{String, Vararg{String}}}="", only_write_changed_cols::Bool=false, 
                          seed::Union{Nothing, Int}=nothing)
    N = lcm(length(adj_case_path), length(adj_case_name), length(adj_case_ext), 
            length(prod_fac), length(add_fac), length(mean), length(sd), 
            length(write_file_path), length(write_file_name))
    adj_case_path = isa(adj_case_path, Union{String, Tuple{String}}) ? match_length(String(adj_case_path), N) : adj_case_path
    adj_case_name = isa(adj_case_name, Union{String, Tuple{String}}) ? match_length(String(adj_case_name), N) : adj_case_name
    adj_case_ext = isa(adj_case_ext, Union{String, Tuple{String}}) ? match_length(String(adj_case_ext), N) : adj_case_ext
    write_file_path = (!overwrite_file & isa(write_file_path, Union{String, Tuple{String}})) ? match_length(String(write_file_path), N) : write_file_path
    write_file_name = (!overwrite_file & isa(write_file_name, Union{String, Tuple{String}})) ? match_length(String(write_file_name), N) : write_file_name

    p_fac = isa(prod_fac, Union{Real, Tuple{Real}}) ? match_length(prod_fac..., N) : prod_fac
    a_fac = isa(add_fac, Union{Real, Tuple{Real}}) ? match_length(add_fac..., N) : add_fac
    m = isa(mean, Union{Real, Tuple{Real}}) ? match_length(mean..., N) : mean
    s = isa(sd, Union{Real, Tuple{Real}}) ? match_length(sd..., N) : sd

    rng = isnothing(seed) ? MersenneTwister() : MersenneTwister(seed)
    sub_seeds = isnothing(seed) ? match_length(seed, N) : Tuple(rand(rng, Int, N))

    @assert length(adj_case_path) == length(adj_case_name) == length(adj_case_ext) == 
            length(write_file_path) == length(write_file_name) ==
            length(p_fac) == length(a_fac) == length(m) == length(s) == length(sub_seeds)

    for idx in 1:N
        adj_params(adj_case_path[idx], adj_case_name[idx], adj_case_ext[idx], c2, c1, c0, p_fac[idx], a_fac[idx], 
        start_x_idx=start_x_idx, end_x_idx=end_x_idx, T=T, mean=m[idx], sd=s[idx], 
        overwrite_file=overwrite_file, write_file_path=write_file_path[idx], write_file_name=write_file_name[idx], 
        only_write_changed_cols=only_write_changed_cols, seed=sub_seeds[idx])
    end
end

# 4
function adj_multi_params(adj_case_path::Union{String, Tuple{String, Vararg{String}}}, adj_case_name::Union{String, Tuple{String, Vararg{String}}}, adj_case_ext::Union{String, Tuple{String, Vararg{String}}}, 
                          c2::Bool, c1::Bool, c0::Bool, vals::Union{VecOrMat{<:Real}, NTuple{N, VecOrMat{<:Real}}}; 
                          start_x_idx::Int=1, end_x_idx::Int=0, T::Type=Float64, mean::Union{Nothing, Real, Tuple{Real, Vararg{Real}}}=nothing, sd::Union{Nothing, Real, Tuple{Real, Vararg{Real}}}=nothing, 
                          overwrite_file::Bool=false, write_file_path::Union{String, Tuple{String, Vararg{String}}}="", write_file_name::Union{String, Tuple{String, Vararg{String}}}="", only_write_changed_cols::Bool=false, 
                          seed::Union{Nothing, Int}=nothing) where {N}
    N = lcm(length(adj_case_path), length(adj_case_name), length(adj_case_ext), 
            length(vals), length(mean), length(sd), 
            length(write_file_path), length(write_file_name))
    adj_case_path = isa(adj_case_path, Union{String, Tuple{String}}) ? match_length(String(adj_case_path), N) : adj_case_path
    adj_case_name = isa(adj_case_name, Union{String, Tuple{String}}) ? match_length(String(adj_case_name), N) : adj_case_name
    adj_case_ext = isa(adj_case_ext, Union{String, Tuple{String}}) ? match_length(String(adj_case_ext), N) : adj_case_ext
    write_file_path = (!overwrite_file & isa(write_file_path, Union{String, Tuple{String}})) ? match_length(String(write_file_path), N) : write_file_path
    write_file_name = (!overwrite_file & isa(write_file_name, Union{String, Tuple{String}})) ? match_length(String(write_file_name), N) : write_file_name

    val = isa(vals, Union{VecOrMat{<:Real}, Tuple{VecOrMat{<:Real}}}) ? match_length(vals..., N) : vals
    m = isa(mean, Union{Real, Tuple{Real}}) ? match_length(mean..., N) : mean
    s = isa(sd, Union{Real, Tuple{Real}}) ? match_length(sd..., N) : sd

    rng = isnothing(seed) ? MersenneTwister() : MersenneTwister(seed)
    sub_seeds = isnothing(seed) ? match_length(seed, N) : Tuple(rand(rng, Int, N))

    @assert length(adj_case_path) == length(adj_case_name) == length(adj_case_ext) == 
            length(write_file_path) == length(write_file_name) ==
            length(val) == length(m) == length(s) == length(sub_seeds)

    for idx in 1:N
        adj_params(adj_case_path[idx], adj_case_name[idx], adj_case_ext[idx], c2, c1, c0, val[idx], 
        start_x_idx=start_x_idx, end_x_idx=end_x_idx, T=T, mean=m[idx], sd=s[idx], 
        overwrite_file=overwrite_file, write_file_path=write_file_path[idx], write_file_name=write_file_name[idx], 
        only_write_changed_cols=only_write_changed_cols, seed=sub_seeds[idx])
    end
end

# 5
function adj_multi_params(adj_case_path::Union{String, Tuple{String, Vararg{String}}}, adj_case_name::Union{String, Tuple{String, Vararg{String}}}, adj_case_ext::Union{String, Tuple{String, Vararg{String}}}, 
                          rateA::Bool, prod_fac::Union{Real, Tuple{Real, Vararg{Real}}}, add_fac::Union{Real, Tuple{Real, Vararg{Real}}}; 
                          start_x_idx::Int=1, end_x_idx::Int=0, T::Type=Float64, mean::Union{Nothing, Real, Tuple{Real, Vararg{Real}}}=nothing, sd::Union{Nothing, Real, Tuple{Real, Vararg{Real}}}=nothing, 
                          overwrite_file::Bool=false, write_file_path::Union{String, Tuple{String, Vararg{String}}}="", write_file_name::Union{String, Tuple{String, Vararg{String}}}="", only_write_changed_cols::Bool=false, 
                          seed::Union{Nothing, Int}=nothing)
    N = lcm(length(adj_case_path), length(adj_case_name), length(adj_case_ext), 
            length(prod_fac), length(add_fac), length(mean), length(sd), 
            length(write_file_path), length(write_file_name))
    adj_case_path = isa(adj_case_path, Union{String, Tuple{String}}) ? match_length(String(adj_case_path), N) : adj_case_path
    adj_case_name = isa(adj_case_name, Union{String, Tuple{String}}) ? match_length(String(adj_case_name), N) : adj_case_name
    adj_case_ext = isa(adj_case_ext, Union{String, Tuple{String}}) ? match_length(String(adj_case_ext), N) : adj_case_ext
    write_file_path = (!overwrite_file & isa(write_file_path, Union{String, Tuple{String}})) ? match_length(String(write_file_path), N) : write_file_path
    write_file_name = (!overwrite_file & isa(write_file_name, Union{String, Tuple{String}})) ? match_length(String(write_file_name), N) : write_file_name

    p_fac = isa(prod_fac, Union{Real, Tuple{Real}}) ? match_length(prod_fac..., N) : prod_fac
    a_fac = isa(add_fac, Union{Real, Tuple{Real}}) ? match_length(add_fac..., N) : add_fac
    m = isa(mean, Union{Real, Tuple{Real}}) ? match_length(mean..., N) : mean
    s = isa(sd, Union{Real, Tuple{Real}}) ? match_length(sd..., N) : sd

    rng = isnothing(seed) ? MersenneTwister() : MersenneTwister(seed)
    sub_seeds = isnothing(seed) ? match_length(seed, N) : Tuple(rand(rng, Int, N))

    @assert length(adj_case_path) == length(adj_case_name) == length(adj_case_ext) == 
            length(write_file_path) == length(write_file_name) ==
            length(p_fac) == length(a_fac) == length(m) == length(s) == length(sub_seeds)

    for idx in 1:N
        adj_params(adj_case_path[idx], adj_case_name[idx], adj_case_ext[idx], rateA, p_fac[idx], a_fac[idx], 
        start_x_idx=start_x_idx, end_x_idx=end_x_idx, T=T, mean=m[idx], sd=s[idx], 
        overwrite_file=overwrite_file, write_file_path=write_file_path[idx], write_file_name=write_file_name[idx], 
        only_write_changed_cols=only_write_changed_cols, seed=sub_seeds[idx])
    end
end

# 6
function adj_multi_params(adj_case_path::Union{String, Tuple{String, Vararg{String}}}, adj_case_name::Union{String, Tuple{String, Vararg{String}}}, adj_case_ext::Union{String, Tuple{String, Vararg{String}}}, 
                          rateA::Bool, vals::Union{VecOrMat{<:Real}, NTuple{N, VecOrMat{<:Real}}}; 
                          start_x_idx::Int=1, end_x_idx::Int=0, T::Type=Float64, mean::Union{Nothing, Real, Tuple{Real, Vararg{Real}}}=nothing, sd::Union{Nothing, Real, Tuple{Real, Vararg{Real}}}=nothing, 
                          overwrite_file::Bool=false, write_file_path::Union{String, Tuple{String, Vararg{String}}}="", write_file_name::Union{String, Tuple{String, Vararg{String}}}="", only_write_changed_cols::Bool=false, 
                          seed::Union{Nothing, Int}=nothing) where {N}
    N = lcm(length(adj_case_path), length(adj_case_name), length(adj_case_ext), 
            length(vals), length(mean), length(sd), 
            length(write_file_path), length(write_file_name))
    adj_case_path = isa(adj_case_path, Union{String, Tuple{String}}) ? match_length(String(adj_case_path), N) : adj_case_path
    adj_case_name = isa(adj_case_name, Union{String, Tuple{String}}) ? match_length(String(adj_case_name), N) : adj_case_name
    adj_case_ext = isa(adj_case_ext, Union{String, Tuple{String}}) ? match_length(String(adj_case_ext), N) : adj_case_ext
    write_file_path = (!overwrite_file & isa(write_file_path, Union{String, Tuple{String}})) ? match_length(String(write_file_path), N) : write_file_path
    write_file_name = (!overwrite_file & isa(write_file_name, Union{String, Tuple{String}})) ? match_length(String(write_file_name), N) : write_file_name

    val = isa(vals, Union{VecOrMat{<:Real}, Tuple{VecOrMat{<:Real}}}) ? match_length(vals..., N) : vals
    m = isa(mean, Union{Real, Tuple{Real}}) ? match_length(mean..., N) : mean
    s = isa(sd, Union{Real, Tuple{Real}}) ? match_length(sd..., N) : sd

    rng = isnothing(seed) ? MersenneTwister() : MersenneTwister(seed)
    sub_seeds = isnothing(seed) ? match_length(seed, N) : Tuple(rand(rng, Int, N))

    @assert length(adj_case_path) == length(adj_case_name) == length(adj_case_ext) == 
            length(write_file_path) == length(write_file_name) ==
            length(val) == length(m) == length(s) == length(sub_seeds)

    for idx in 1:N
        adj_params(adj_case_path[idx], adj_case_name[idx], adj_case_ext[idx], rateA, val[idx], 
        start_x_idx=start_x_idx, end_x_idx=end_x_idx, T=T, mean=m[idx], sd=s[idx], 
        overwrite_file=overwrite_file, write_file_path=write_file_path[idx], write_file_name=write_file_name[idx], 
        only_write_changed_cols=only_write_changed_cols, seed=sub_seeds[idx])
    end
end

# 7
function adj_multi_params(adj_case_path::Union{String, Tuple{String, Vararg{String}}}, adj_case_name::Union{String, Tuple{String, Vararg{String}}}, adj_case_ext::Tuple{String, Vararg{String}}, 
                          P::Union{Bool, Tuple{Bool, Vararg{Bool}}}, Q::Union{Bool, Tuple{Bool, Vararg{Bool}}}, c2::Bool, c1::Bool, c0::Bool, rateA::Bool, prod_fac::Union{Real, Tuple{Real, Vararg{Real}}}, add_fac::Union{Real, Tuple{Real, Vararg{Real}}}; 
                          start_x_idx::Int=1, end_x_idx::Int=0, T::Type=Float64, mean::Union{Nothing, Real, Tuple{Real, Vararg{Real}}}=nothing, sd::Union{Nothing, Real, Tuple{Real, Vararg{Real}}}=nothing, 
                          overwrite_file::Bool=false, write_file_path::Union{String, Tuple{String, Vararg{String}}}="", write_file_name::Union{String, Tuple{String, Vararg{String}}}="", only_write_changed_cols::Bool=false, 
                          seed::Union{Nothing, Int}=nothing)
    N = lcm(length(adj_case_path), length(adj_case_name), length(adj_case_ext), 
            length(prod_fac), length(add_fac), length(mean), length(sd), 
            length(write_file_path), length(write_file_name))
    adj_case_path = isa(adj_case_path, Union{String, Tuple{String}}) ? match_length(String(adj_case_path), N) : adj_case_path
    adj_case_name = isa(adj_case_name, Union{String, Tuple{String}}) ? match_length(String(adj_case_name), N) : adj_case_name
    write_file_path = (!overwrite_file & isa(write_file_path, Union{String, Tuple{String}})) ? match_length(String(write_file_path), N) : write_file_path
    write_file_name = (!overwrite_file & isa(write_file_name, Union{String, Tuple{String}})) ? match_length(String(write_file_name), N) : write_file_name

    p_fac = isa(prod_fac, Union{Real, Tuple{Real}}) ? match_length(prod_fac..., N) : prod_fac
    a_fac = isa(add_fac, Union{Real, Tuple{Real}}) ? match_length(add_fac..., N) : add_fac
    m = isa(mean, Union{Real, Tuple{Real}}) ? match_length(mean..., N) : mean
    s = isa(sd, Union{Real, Tuple{Real}}) ? match_length(sd..., N) : sd

    rng = isnothing(seed) ? MersenneTwister() : MersenneTwister(seed)
    sub_seeds = isnothing(seed) ? match_length(seed, N) : Tuple(rand(rng, Int, N))

    @assert length(adj_case_path) == length(adj_case_name) == length(adj_case_ext) == 
            length(write_file_path) == length(write_file_name) ==
            length(p_fac) == length(a_fac) == length(m) == length(s) == length(sub_seeds)

    for idx in 1:N
        adj_params(adj_case_path[idx], adj_case_name[idx], adj_case_ext[idx], P, Q, c2, c1, c0, rateA, p_fac[idx], a_fac[idx], 
        start_x_idx=start_x_idx, end_x_idx=end_x_idx, T=T, mean=m[idx], sd=s[idx], 
        overwrite_file=overwrite_file, write_file_path=write_file_path[idx], write_file_name=write_file_name[idx], 
        only_write_changed_cols=only_write_changed_cols, seed=sub_seeds[idx])
    end
end

# 8
function adj_multi_params(adj_case_path::Union{String, Tuple{String, Vararg{String}}}, adj_case_name::Union{String, Tuple{String, Vararg{String}}}, adj_case_ext::NTuple{N, String}, 
                          P::Union{Bool, Tuple{Bool, Vararg{Bool}}}, Q::Union{Bool, Tuple{Bool, Vararg{Bool}}}, c2::Bool, c1::Bool, c0::Bool, rateA::Bool, vals::Union{VecOrMat{<:Real}, NTuple{N, VecOrMat{<:Real}}}; 
                          start_x_idx::Int=1, end_x_idx::Int=0, T::Type=Float64, mean::Union{Nothing, Real, Tuple{Real, Vararg{Real}}}=nothing, sd::Union{Nothing, Real, Tuple{Real, Vararg{Real}}}=nothing, 
                          overwrite_file::Bool=false, write_file_path::Union{String, Tuple{String, Vararg{String}}}="", write_file_name::Union{String, Tuple{String, Vararg{String}}}="", only_write_changed_cols::Bool=false, 
                          seed::Union{Nothing, Int}=nothing) where {N}
    N = lcm(length(adj_case_path), length(adj_case_name), length(adj_case_ext), 
            length(vals), length(mean), length(sd), 
            length(write_file_path), length(write_file_name))
    adj_case_path = isa(adj_case_path, Union{String, Tuple{String}}) ? match_length(String(adj_case_path), N) : adj_case_path
    adj_case_name = isa(adj_case_name, Union{String, Tuple{String}}) ? match_length(String(adj_case_name), N) : adj_case_name
    write_file_path = (!overwrite_file & isa(write_file_path, Union{String, NTuple{1, String}})) ? match_length(String(write_file_path), N) : write_file_path
    write_file_name = (!overwrite_file & isa(write_file_name, Union{String, Tuple{String}})) ? match_length(String(write_file_name), N) : write_file_name

    val = isa(vals, Union{VecOrMat{<:Real}, Tuple{VecOrMat{<:Real}}}) ? match_length(vals..., N) : vals
    m = isa(mean, Union{Real, Tuple{Real}}) ? match_length(mean..., N) : mean
    s = isa(sd, Union{Real, Tuple{Real}}) ? match_length(sd..., N) : sd

    rng = isnothing(seed) ? MersenneTwister() : MersenneTwister(seed)
    sub_seeds = isnothing(seed) ? match_length(seed, N) : Tuple(rand(rng, Int, N))

    @assert length(adj_case_path) == length(adj_case_name) == length(adj_case_ext) == 
            length(write_file_path) == length(write_file_name) ==
            length(val) == length(m) == length(s) == length(sub_seeds)

    for idx in 1:N
        adj_params(adj_case_path[idx], adj_case_name[idx], adj_case_ext[idx], P, Q, c2, c1, c0, rateA, val[idx], 
        start_x_idx=start_x_idx, end_x_idx=end_x_idx, T=T, mean=m[idx], sd=s[idx], 
        overwrite_file=overwrite_file, write_file_path=write_file_path[idx], write_file_name=write_file_name[idx], 
        only_write_changed_cols=only_write_changed_cols, seed=sub_seeds[idx])
    end
end
