# Balance the slack of Pg from xbar after perturbing Pg

function get_xbar_dict(cascades_root::String, case_name::String, oppt_dir_name::String)
    path_dict = get_path_dict(cascades_root, case_name)
    xbar_dict = Dict{Symbol, Array}()
    for file in readdir(path_dict[:operatingdata_path] * oppt_dir_name * "/")
        file_name = replace(file, ".csv" => "")
        arr = readdlm(path_dict[:operatingdata_path] * oppt_dir_name * "/" * file)
        xbar_dict[Symbol(file_name)] = size(arr, 2) == 1 ? vec(arr) : arr
    end
    return xbar_dict
end

function get_xbar_arr(cascades_root::String, case_name::String, oppt_dir_name::String, sol_file_name::String)
    xbar_dict = get_xbar_dict(cascades_root, case_name, oppt_dir_name)
    xbar_arr = xbar_dict[Symbol(sol_file_name)]
    return xbar_arr
end

function get_perturbed_xbar_arr(cascades_root::String, case_name::String, oppt_dir_name::String, sol_file_name::String; 
                                start_x_idx::Int=1, y_idx::Union{Nothing, Vector{Int}}=nothing, 
                                mean::Union{Nothing, Int, Float64}=nothing, sd::Union{Nothing, Int, Float64}=nothing, 
                                scale_fac::Union{Nothing, Int, Float64}=nothing,
                                discard_neg_vals::Bool=true, seed::Union{Nothing, Int}=nothing)
    if !isnothing(scale_fac)
        @assert scale_fac >= 0
    end

    xbar_arr = get_xbar_arr(cascades_root, case_name, oppt_dir_name, sol_file_name)
    if (isnothing(mean) & isnothing(sd)) & isnothing(scale_fac)
        perturbed_xbar_arr = xbar_arr 
    elseif (!isnothing(mean) & !isnothing(sd)) & isnothing(scale_fac)
        perturbed_xbar_arr = add_gaussian_noise(xbar_arr, mean, sd, seed)
    elseif (isnothing(mean) & isnothing(sd)) & !isnothing(scale_fac)
        perturbed_xbar_arr = xbar_arr .* scale_fac
    else 
        throw(DomainError("Perturbation parameters are not properly defined."))
    end 

    abs_perturbed_xbar_arr = discard_neg_vals ? undo_neg_vals(perturbed_xbar_arr, start_x_idx, y_idx, xbar_arr) : perturbed_xbar_arr
    return abs_perturbed_xbar_arr
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

function balance_slack(cascades_root::String, case_name::String, oppt_dir_name::String, sol_file_name::String, file_name::String;
                       start_x_idx::Int=1, y_idx::Union{Nothing, Vector{Int}}=nothing, 
                       mean::Union{Nothing, Int, Float64}=nothing, sd::Union{Nothing, Int, Float64}=nothing, 
                       scale_fac::Union{Nothing, Int, Float64}=nothing,
                       discard_neg_vals::Bool=true, seed::Union{Nothing, Int}=nothing)
    slack_bus_type = 3
    slack_id = Int(first(get_bus_id(cascades_root, case_name, file_name, slack_bus_type)))
    gen_arr = get_case_arr(cascades_root, case_name, file_name, ".gen")
    slack_idx = findall(x -> x == slack_id, gen_arr[:, 1])
    @assert length(slack_idx) == 1
    slack_idx = first(slack_idx)

    xbar_arr = get_xbar_arr(cascades_root, case_name, oppt_dir_name, sol_file_name)
    perturbed_xbar_arr = get_perturbed_xbar_arr(cascades_root, case_name, oppt_dir_name, sol_file_name,
                                                start_x_idx=start_x_idx, y_idx=y_idx, mean=mean, sd=sd, 
                                                scale_fac=scale_fac, discard_neg_vals=discard_neg_vals, seed=seed)

    perturbed_slack = perturbed_xbar_arr[slack_idx]
    perturbed_sum = sum(perturbed_xbar_arr - xbar_arr)
    if perturbed_slack < perturbed_sum
        perturbed_nonslack_sum = perturbed_sum - (perturbed_slack - xbar_arr[slack_idx])
        scale_fac = perturbed_slack / perturbed_nonslack_sum
        perturbed_xbar_arr .*= scale_fac
        perturbed_sum = sum(perturbed_xbar_arr - xbar_arr)
    end

    perturbed_xbar_arr[slack_idx] -= perturbed_sum
    @assert sum(perturbed_xbar_arr) - sum(xbar_arr) <= eps(Float32)
    return perturbed_xbar_arr
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
