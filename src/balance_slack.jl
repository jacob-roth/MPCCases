# Balance the slack of Pg from xbar after perturbing Pg

function get_perturbed_xbar_arr(cascades_root::String, case_name::String, oppt_dir_name::String, sol_file_name::String; 
                                start_x_idx::Int=1, y_idx::Union{Nothing, Vector{Int}}=nothing, 
                                mean::Union{Nothing, Int, Float64}=nothing, sd::Union{Nothing, Int, Float64}=nothing, 
                                scale_fac::Union{Nothing, Int, Float64}=nothing, max_delta_prop::Union{Nothing, Int, Float64}=nothing,
                                discard_neg_vals::Bool=true, idx_arr::Union{Nothing, Vector{Int}}=nothing,
                                seed::Union{Nothing, Int}=nothing)
    if !isnothing(scale_fac)
        @assert scale_fac >= 0
    end

    xbar_arr = get_xbar_arr(cascades_root, case_name, oppt_dir_name, sol_file_name)

    if (isnothing(mean) & isnothing(sd)) & isnothing(scale_fac) & isnothing(max_delta_prop)
        perturbed_xbar_arr = xbar_arr 
    elseif (!isnothing(mean) & !isnothing(sd)) & isnothing(scale_fac) & isnothing(max_delta_prop)
        perturbed_xbar_arr = add_gaussian_noise(xbar_arr, mean, sd, seed)
    elseif (isnothing(mean) & isnothing(sd)) & !isnothing(scale_fac) & isnothing(max_delta_prop)
        perturbed_xbar_arr = xbar_arr .* scale_fac
    elseif (isnothing(mean) & isnothing(sd)) & !isnothing(scale_fac) & !isnothing(max_delta_prop)
        agg_xbar_arr = sum_xbar_arr(xbar_arr, idx_arr=idx_arr)
        max_delta_val = max_delta_prop * agg_xbar_arr
        scaled_xbar_arr = xbar_arr .* scale_fac
        perturbed_xbar_arr = scaled_xbar_arr .* (scaled_xbar_arr .<= max_delta_val) + (max_delta_val * (scaled_xbar_arr .>= max_delta_val))
    elseif (isnothing(mean) & isnothing(sd)) & isnothing(scale_fac) & !isnothing(max_delta_prop)
        agg_xbar_arr = sum_xbar_arr(xbar_arr, idx_arr=idx_arr)
        max_delta_val = max_delta_prop * agg_xbar_arr
        perturbed_xbar_arr = xbar_arr .* (xbar_arr .<= max_delta_val) + (max_delta_val * (xbar_arr .>= max_delta_val))
    else 
        throw(DomainError("Perturbation parameters are not properly defined."))
    end 

    abs_perturbed_xbar_arr = discard_neg_vals ? undo_neg_vals(perturbed_xbar_arr, start_x_idx, y_idx, xbar_arr) : perturbed_xbar_arr
    return abs_perturbed_xbar_arr
end

function sum_xbar_arr(xbar_arr::VecOrMat{<:Real}; idx_arr::Union{Nothing, Vector{Int}}=nothing)
    if isa(xbar_arr, Vector{<:Real})
        return sum(xbar_arr)
    elseif !isnothing(idx_arr)
        return sum(xbar_arr[idx_arr])
    else
        throw(DomainError("idx_arr is not properly defined."))
    end
end

function balance_slack(cascades_root::String, case_name::String, oppt_dir_name::String, sol_file_name::String, file_name::String;
                       start_x_idx::Int=1, y_idx::Union{Nothing, Vector{Int}}=nothing, 
                       mean::Union{Nothing, Int, Float64}=nothing, sd::Union{Nothing, Int, Float64}=nothing, 
                       scale_fac::Union{Nothing, Int, Float64}=nothing, max_delta_prop::Union{Nothing, Int, Float64}=nothing,
                       discard_neg_vals::Bool=true, idx_arr::Union{Nothing, Vector{Int}}=nothing, 
                       seed::Union{Nothing, Int}=nothing, err::Real=1e-9)
    slack_bus_type = 3
    slack_id = Int(first(get_bus_id(cascades_root, case_name, file_name, slack_bus_type)))
    gen_arr = get_case_arr(cascades_root, case_name, file_name, ".gen")
    slack_idx = findall(x -> x == slack_id, gen_arr[:, 1])
    @assert length(slack_idx) == 1
    slack_idx = first(slack_idx)

    xbar_arr = get_xbar_arr(cascades_root, case_name, oppt_dir_name, sol_file_name)
    perturbed_xbar_arr = get_perturbed_xbar_arr(cascades_root, case_name, oppt_dir_name, sol_file_name,
                                                start_x_idx=start_x_idx, y_idx=y_idx, mean=mean, sd=sd, 
                                                scale_fac=scale_fac, max_delta_prop=max_delta_prop,
                                                discard_neg_vals=discard_neg_vals, idx_arr=idx_arr,
                                                seed=seed)

    perturbed_slack = perturbed_xbar_arr[slack_idx]
    perturbed_sum = sum(perturbed_xbar_arr - xbar_arr)
    if perturbed_slack < perturbed_sum
        perturbed_nonslack_sum = perturbed_sum - (perturbed_slack - xbar_arr[slack_idx])
        scale_fac = perturbed_slack / perturbed_nonslack_sum
        perturbed_xbar_arr .*= scale_fac
        perturbed_sum = sum(perturbed_xbar_arr - xbar_arr)
    end

    perturbed_xbar_arr[slack_idx] -= perturbed_sum
    @assert abs(sum(perturbed_xbar_arr) - sum(xbar_arr)) <= err
    return perturbed_xbar_arr
end

function get_bus_id(cascades_root::String, case_name::String, file_name::String, bus_type::Int)
    path_dict = get_path_dict(cascades_root, case_name)
    bus_arr = readdlm(path_dict[:casedata_path] * file_name * ".bus")
    bus_idx = findall(x -> x == bus_type, bus_arr[:, 2])
    if bus_type == 3
        @assert length(bus_idx) == 1
    end
    bus_id = bus_arr[bus_idx, 1]
    return bus_id
end

function write_perturbed_Pg_file(cascades_root::String, case_name::String, oppt_dir_name::String, 
                                 write_file_name::String, perturbed_xbar_arr::VecOrMat{<:Real})
    path_dict = get_path_dict(cascades_root, case_name)
    operatingdata_path = path_dict[:operatingdata_path]
    open(operatingdata_path * oppt_dir_name * "/" * write_file_name * ".gen", "w") do io
        writedlm(io, perturbed_xbar_arr)
    end
end

function get_true_scale_fac(perturbed_xbar_arr::VecOrMat{<:Real}, xbar_arr::VecOrMat{<:Real})
    @assert size(perturbed_xbar_arr) == size(xbar_arr)
    scale_fac = perturbed_xbar_arr ./ xbar_arr
    return scale_fac
end
