# Perturb Pg from xbar

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
                                discard_neg_vals::Bool=true, seed::Union{Nothing, Int}=nothing)
    xbar_arr = get_xbar_arr(cascades_root, case_name, oppt_dir_name, sol_file_name)
    perturbed_xbar_arr = isnothing(mean) | isnothing(sd) ? xbar_arr : add_gaussian_noise(xbar_arr, mean, sd, seed)
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
    bus_id = bus_arr[1, bus_idx]
    return bus_id
end

function balance_slack(cascades_root::String, case_name::String, oppt_dir_name::String, sol_file_name::String, file_name::String;
                       start_x_idx::Int=1, y_idx::Union{Nothing, Vector{Int}}=nothing, 
                       mean::Union{Nothing, Int, Float64}=nothing, sd::Union{Nothing, Int, Float64}=nothing, 
                       discard_neg_vals::Bool=true, seed::Union{Nothing, Int}=nothing)
    slack_bus_type = 3
    slack_id = get_bus_id(cascades_root, case_name, file_name, slack_bus_type)
    gen_arr = get_case_arr(cascades_root, case_name, file_name, ".gen")
    slack_idx = findall(x -> x == slack_id, gen_arr[:, 1])
    @assert length(slack_idx) == 1

    xbar_arr = get_xbar_arr(cascades_root, case_name, oppt_dir_name, sol_file_name)
    perturbed_xbar_arr = get_perturbed_xbar_arr(cascades_root, case_name, oppt_dir_name, sol_file_name,
                                                start_x_idx=start_x_idx, y_idx=y_idx, mean=mean, sd=sd, 
                                                discard_neg_vals=discard_neg_vals, seed=seed)
    perturbation_sum = sum(perturbed_xbar_arr - xbar_arr)
    perturbed_xbar_arr[slack_idx, 1] -= perturbation_sum
    
end

