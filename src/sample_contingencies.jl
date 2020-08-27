function get_num_lines(casedata::CaseData)
    lines = casedata.opf.lines
    @assert lines == unique(lines)
    return length(lines)
end

function get_num_buses(casedata::CaseData)
    buses = casedata.opf.buses
    @assert buses == unique(buses)
    return length(buses)
end

function get_unif_rv(a::Union{Nothing, Int}=nothing, b::Union{Nothing, Int}=nothing; 
                     dims::Union{Nothing, Int, NTuple{N, Int}}=nothing, 
                     seed::Union{Nothing, Int}=nothing) where {N}
    rng = isnothing(seed) ? MersenneTwister() : MersenneTwister(seed)
    unif_rv =   isnothing(a) & isnothing(b) ? rand(rng) :
                isnothing(dims) ? rand(rng, a:b) : 
                rand(rng, a:b, dims)
    return unif_rv
end

function get_initial_failed_line_id(casedata::CaseData; 
                                    dims::Union{Nothing, Int, NTuple{N, Int}}=nothing, 
                                    seed::Union{Nothing, Int}=nothing) where {N}
    num_lines = get_num_lines(casedata)
    failed_line_id = get_unif_rv(1, num_lines, dims=dims, seed=seed)
    return failed_line_id
end

function get_zipf_pmf(s::Real, N::Int, k::Union{UnitRange{Int}, Int})
    @assert (s >= 0) & (N > 0)
    if isa(k, Int)
        @assert (1 <= k <= N)
        zipf_pmf = (1 / (k^s)) / sum([1/(n^s) for n in 1:N])
    else
        [@assert (1 <= l <= N) for l in k]
        denom = sum([1/(n^s) for n in 1:N])
        zipf_pmf = [(1 / (l^s)) / denom for l in k]
    end
    return zipf_pmf
end

function get_zipf_cdf(s::Real, N::Int, k::Union{UnitRange{Int}, Int})
    if isa(k, Int)
        zipf_cdf = sum([get_zipf_pmf(s, N, l) for l in 1:k])
    else
        zipf_cdf = [sum([get_zipf_pmf(s, N, m) for m in 1:l]) for l in k]
    end
    return zipf_cdf
end

function get_zipf_dict(s::Real, N::Int, k::Union{UnitRange{Int}, Int}; pmf_flag::Bool=true)
    zipf_arr = pmf_flag ? get_zipf_pmf(s, N, k) : get_zipf_cdf(s, N, k)
    zipf_dict = Dict{Int, Union{Int,Float64}}()
    if isa(k, Int)
        zipf_dict[k] = zipf_arr
    else
        for l in k
            zipf_dict[l] = zipf_arr[l]
        end
    end
    return zipf_dict
end

function get_bus_ids(casedata::CaseData, line_id::Int)
    lines = casedata.opf.lines
    from_bus = first(lines[lines.id .== line_id].from)
    to_bus = first(lines[lines.id .== line_id].to)
    if from_bus < to_bus
        return (from_bus, to_bus)
    else
        return (to_bus, from_bus)
    end
end

function get_line_bus_pairs(casedata::CaseData; reverse_dict::Bool=false)
    lines = casedata.opf.lines
    line_bus_pairs = Dict{Int, Tuple{Int, Int}}()
    for line_id in lines.id
        line_bus_pairs[line_id] = get_bus_ids(casedata, line_id)
    end
    if reverse_dict
        return Dict(value => key for (key, value) in line_bus_pairs)
    else
        return line_bus_pairs
    end
end

function get_line_bus_pairs(casedata::CaseData, line_id::Union{Int, Array{Int}}; reverse_dict::Bool=false)
    line_bus_pairs = Dict{Int, Tuple{Int, Int}}()
    for l_id in line_id
        line_bus_pairs[l_id] = get_bus_ids(casedata, l_id)
    end
    if reverse_dict
        return Dict(value => key for (key, value) in line_bus_pairs)
    else
        return line_bus_pairs
    end
end

function add_bus_to_tree(bus_tree::Array{Union{Nothing, Array{Int}}}, from_bus::Int, to_bus::Int)
    if isnothing(bus_tree[from_bus])
        return [to_bus]
    else
        return cat(bus_tree[from_bus], to_bus, dims=1)
    end
end

function add_bus_as_neighbour(neigbour_dict::Dict{Int, Array{Tuple{Int, Int}}}, distance_key::Int, target_bus_id::Int, parent_bus_id::Int)
    bus_pair = parent_bus_id < target_bus_id ? (parent_bus_id, target_bus_id) : (target_bus_id, parent_bus_id)
    if !(haskey(neigbour_dict, distance_key))
        return [bus_pair]
    else
        return unique(cat(neigbour_dict[distance_key], bus_pair, dims=1))
    end
end

function get_bus_tree(casedata::CaseData)
    num_buses = get_num_buses(casedata)
    bus_pairs = get_line_bus_pairs(casedata)
    bus_tree = Array{Union{Nothing, Array{Int}}}(nothing, (num_buses,1))
    for bus_pair in values(bus_pairs)
        from_bus, to_bus = bus_pair
        bus_tree[from_bus] = add_bus_to_tree(bus_tree, from_bus, to_bus)
        bus_tree[to_bus] = add_bus_to_tree(bus_tree, to_bus, from_bus)
    end
    return bus_tree
end

function dfs(distance::Int, target_bus_id::Int, bus_tree::Array{Union{Nothing, Array{Int}}};
             neighbour_dict::Dict{Int, Array{Tuple{Int, Int}}}=Dict{Int, Array{Tuple{Int, Int}}}(), 
             parent_bus_id::Int=-1, depth::Int=0)
    if distance <= -1
        return nothing
    end
    
    neighbour_dict[depth] = add_bus_as_neighbour(neighbour_dict, depth, target_bus_id, parent_bus_id)
    for child_bus_id in bus_tree[target_bus_id]
        if child_bus_id != parent_bus_id
            depth += 1
            dfs(distance - 1, child_bus_id, bus_tree, 
                neighbour_dict=neighbour_dict, parent_bus_id=target_bus_id, depth=depth)
            depth -= 1
        end
    end
    return neighbour_dict
end

function remove_cycled_lines(neighbour_dict::Dict{Int, Array{Tuple{Int, Int}}})
    neighbour_dict_cp = copy(neighbour_dict)
    tmp_arr = []
    for distance in sort(collect(keys(neighbour_dict_cp)))
        for child_bus in neighbour_dict_cp[distance]
            if child_bus ∉ tmp_arr
                tmp_arr = cat(tmp_arr, child_bus, dims=1)
            else
                neighbour_dict_cp[distance] = filter(x -> x ≠ child_bus, neighbour_dict_cp[distance])
            end
        end
    end
    return neighbour_dict_cp
end

function remove_cycled_lines(neighbour_dict::Dict{Int, Array{Int}})
    neighbour_dict_cp = copy(neighbour_dict)
    tmp_arr = []
    for distance in sort(collect(keys(neighbour_dict_cp)))
        for child_bus in neighbour_dict_cp[distance]
            if child_bus ∉ tmp_arr
                tmp_arr = cat(tmp_arr, child_bus, dims=1)
            else
                neighbour_dict_cp[distance] = filter(x -> x ≠ child_bus, neighbour_dict_cp[distance])
            end
        end
    end
    return neighbour_dict_cp
end

function get_children_of_failed_line(casedata::CaseData, failed_line_id::Int, distance::Int;
                                     recursive::Bool=true)
    failed_line = get_line_bus_pairs(casedata, failed_line_id, reverse_dict=false)
    bus_tree = get_bus_tree(casedata)

    children_dict = Dict{Tuple{Int, Int, Int}, Array{Tuple{Int,Int}}}()
    for failed_bus in collect(Iterators.flatten(values(failed_line)))
        neighbour_dict = remove_cycled_lines(dfs(distance, failed_bus, bus_tree))
        if recursive
            for depth in setdiff(keys(neighbour_dict), 0)
                children_dict[(failed_line_id, failed_bus, depth)] = [neighbour_dict[depth][idx] for idx in eachindex(neighbour_dict[depth])]
            end
        else
            children_dict[(failed_line_id, failed_bus, distance)] = [neighbour_dict[distance][idx] for idx in eachindex(neighbour_dict[distance])]
        end
    end
    return children_dict
end

function key_to_tuple(string::String)
    endless_string = filter(x -> x ∉ Set(['(', ')']), string)
    string_arr = [parse(Int, str) for str in split(endless_string, ',')]
    string_tuple = Tuple(string_arr)
    return string_tuple
end

function value_to_tuple_array(value::Array{Any})
    val_arr = [Tuple{Int, Int}(val_pair) for val_pair in value]
    return val_arr
end

function convert_children_dict(json_dict::Dict{String, Any})
    children_dict = Dict{Tuple{Int, Int, Int}, Array{Tuple{Int, Int}}}()
    for (key, value) in json_dict
        key_tuple = key_to_tuple(key)
        value_arr = value_to_tuple_array(value)
        children_dict[key_tuple] = value_arr
    end
    return children_dict
end

function load_children_dict(path_to_children_dict::String)
    json_dict = Dict()
    open(path_to_children_dict, "r") do io
        json_file = read(io, String)
        json_dict = JSON.parse(json_file)
    end
    return json_dict
end

function add_to_loaded_dict(casedata::CaseData, initial_failed_line_id::Int, distance::Int,
                            loaded_children_dict::Dict{Tuple{Int, Int, Int}, Array{Tuple{Int, Int}}})
    loaded_children_dict_cp = copy(loaded_children_dict)
    tmp_dict = get_children_of_failed_line(casedata, initial_failed_line_id, distance, recursive=true)
    for key in keys(tmp_dict)
        if !(haskey(loaded_children_dict_cp, key))
            loaded_children_dict_cp[key] = tmp_dict[key]
        end
    end
    return loaded_children_dict_cp
end

function cat_to_existing_dict(curr_children_dict::Dict{Tuple{Int, Int, Int}, Array{Tuple{Int, Int}}}, 
                              path_to_children_dict::Union{Nothing, String})
    if isnothing(path_to_children_dict)
        return curr_children_dict
    end
    existing_children_dict = convert_children_dict(load_children_dict(path_to_children_dict))
    for key in keys(curr_children_dict)
        if !(haskey(existing_children_dict, key))
            existing_children_dict[key] = curr_children_dict[key]
        end
    end
    return existing_children_dict
end

function save_children_dict(curr_children_dict::Dict{Tuple{Int, Int, Int}, Array{Tuple{Int, Int}}}, 
                            path_to_children_dict::Union{Nothing, String}; 
                            overwrite_file::Bool=true, write_file_path::Union{Nothing, String}=nothing)
    children_dict = cat_to_existing_dict(curr_children_dict, path_to_children_dict)
    json_data = JSON.json(children_dict)
    if overwrite_file
        open(path_to_children_dict, "w") do io
            write(io, json_data)
        end
    else
        if !(isnothing(write_file_path))
            open(write_file_path, "w") do io
                write(io, json_data)
            end
        else
            throw(UndefKeywordError("write_file_path not defined"))
        end
    end
end

function remove_already_failed_line(casedata::CaseData, children_dict::Dict{Tuple{Int, Int, Int}, Array{Tuple{Int, Int}}}, failed_line_id::Int)
    failed_bus_pair = first(values(get_line_bus_pairs(casedata, failed_line_id, reverse_dict=false)))
    for (key, line_arr) in children_dict
        if failed_bus_pair ∈ line_arr
            children_dict[key] = filter(x -> x ≠ failed_bus_pair, line_arr)
        end
    end
    return children_dict
end

function match_to_line_id(children_dict::Dict{Tuple{Int, Int, Int}, Array{Tuple{Int, Int}}}, rev_line_bus_pairs::Dict{Tuple{Int, Int}, Int})
    line_id_dict = Dict{Int, Array{Int}}()
    for (line_id, bus_id, depth) in keys(children_dict)
        for bus_pair in children_dict[(line_id, bus_id, depth)]
            if !(haskey(line_id_dict, depth))
                line_id_dict[depth] = [rev_line_bus_pairs[bus_pair]]
            else
                line_id_dict[depth] = unique(cat(line_id_dict[depth], rev_line_bus_pairs[bus_pair], dims=1))
            end
        end
    end
    return line_id_dict
end

function get_second_failed_line_id(casedata::CaseData, initial_failed_line_id::Int, s::Real, distance::Int;
                                   recursive::Bool=true, seed::Union{Nothing, Int}=nothing, load_dict::Bool=true, save_dict::Bool=false, 
                                   path_to_children_dict::Union{Nothing, String}=nothing, overwrite_file::Bool=true, write_file_path::Union{Nothing, String}=nothing)
    if load_dict
        @assert !(isnothing(path_to_children_dict))
        children_dict = convert_children_dict(load_children_dict(path_to_children_dict))
        bus_ids = get_bus_ids(casedata, initial_failed_line_id)
        for bus_id in bus_ids
            if (initial_failed_line_id, bus_id, distance) ∉ keys(children_dict)
                children_dict = add_to_loaded_dict(casedata, initial_failed_line_id, distance, children_dict)
            end
        end
    else
        children_dict = get_children_of_failed_line(casedata, initial_failed_line_id, distance, recursive=recursive)
    end
    
    if save_dict
        save_children_dict(children_dict, path_to_children_dict, overwrite_file=overwrite_file, write_file_path=write_file_path)
    end

    valid_children_dict = remove_already_failed_line(casedata, children_dict, initial_failed_line_id)
    rev_line_bus_pairs = get_line_bus_pairs(casedata, reverse_dict=true)
    valid_lines = remove_cycled_lines(match_to_line_id(valid_children_dict, rev_line_bus_pairs))

    [@assert (1 <= k <= distance) for k in keys(valid_lines)]
    cdf_prob_dict = get_zipf_dict(s, distance, 1:distance, pmf_flag=false)

    rng = isnothing(seed) ? MersenneTwister() : MersenneTwister(seed)
    distance_seed, line_seed = isnothing(seed) ? (nothing, nothing) : Tuple(abs.(rand(rng, Int, 2)))

    distance_unif_rv = get_unif_rv(seed=distance_seed)
    failing_bus_distance = 0
    for depth in sort(collect(keys(cdf_prob_dict)))
        if distance_unif_rv >= cdf_prob_dict[depth]
            failing_bus_distance = depth
        end
    end
    failing_bus_distance += 1

    candidate_lines = valid_lines[failing_bus_distance]
    num_candidate_lines = length(candidate_lines)
    line_unif_rv = get_unif_rv(1, num_candidate_lines, seed=line_seed)
    second_failed_line_id = candidate_lines[line_unif_rv]
    return second_failed_line_id
end

function get_second_failed_line_id(casedata::CaseData, initial_failed_line_id::Array{Int}, s::Real, distance::Int;
                                   recursive::Bool=true, seed::Union{Nothing, Int}=nothing, load_dict::Bool=true, save_dict::Bool=false, 
                                   path_to_children_dict::Union{Nothing, String}=nothing, overwrite_file::Bool=true, write_file_path::Union{Nothing, String}=nothing)
    second_failed_line_id = similar(initial_failed_line_id)
    num_sets, num_cores = size(initial_failed_line_id, 1), size(initial_failed_line_id, 2)
    for sets_idx in 1:num_sets
        for cores_idx in 1:num_cores
            failed_line_id = initial_failed_line_id[sets_idx, cores_idx]
            second_failed_line_id[sets_idx, cores_idx] = get_second_failed_line_id(casedata, failed_line_id, s, distance, 
                                                                                   recursive=recursive, seed=seed, load_dict=load_dict, save_dict=save_dict,
                                                                                   path_to_children_dict=path_to_children_dict, overwrite_file=overwrite_file, 
                                                                                   write_file_path=write_file_path)
        end
    end
    return second_failed_line_id
end