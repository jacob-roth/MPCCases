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

function get_unif_rv(a::Int, b::Int; 
                     dims::Union{Nothing, Int, NTuple{N, Int}}=nothing, 
                     seed::Union{Nothing, Int}=nothing) where {N}
    rng = isnothing(seed) ? MersenneTwister() : MersenneTwister(seed)
    unif_rv = isnothing(dims) ? rand(rng, a:b) : rand(rng, a:b, dims)
    return unif_rv
end

function get_initial_failed_line_id(casedata::CaseData; 
                                    dims::Union{Nothing, Int, NTuple{N, Int}}=nothing, 
                                    seed::Union{Nothing, Int}=nothing) where {N}
    num_lines = get_num_lines(casedata)
    failed_line_id = get_unif_rv(1, num_lines, dims=dims, seed=seed)
    return failed_line_id
end

function get_zipf_pmf(s::Union{Int,Float64}, N::Int, k::Union{UnitRange{Int},Int})
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

function get_zipf_cdf(s::Union{Int,Float64}, N::Int, k::Union{UnitRange{Int},Int})
    if isa(k, Int)
        zipf_cdf = sum([get_zipf_pmf(s, N, l) for l in 1:k])
    else
        zipf_cdf = [sum([get_zipf_pmf(s, N, m) for m in 1:l]) for l in k]
    end
    return zipf_cdf
end

function add_bus_to_tree(bus_tree::Array{Union{Nothing, Array{Int}}}, from_bus::Int, to_bus::Int)
    if isnothing(bus_tree[from_bus])
        return [to_bus]
    else
        return cat(bus_tree[from_bus], to_bus, dims=1)
    end
end

function add_bus_as_neighbour(neigbour_dict::Dict{Int, Array{Int}}, 
                              distance_key::Int, 
                              target_bus_id::Int)
    if !(haskey(neigbour_dict, distance_key))
        return [target_bus_id]
    else
        return unique(cat(neigbour_dict[distance_key], target_bus_id, dims=1))
    end
end

function get_bus_tree(casedata::CaseData)
    num_buses = get_num_buses(casedata)
    bus_pairs = get_bus_pairs(casedata)
    bus_tree = Array{Union{Nothing, Array{Int}}}(nothing, (num_buses,1))
    for bus_pair in values(bus_pairs)
        from_bus, to_bus = bus_pair
        bus_tree[from_bus] = add_bus_to_tree(bus_tree, from_bus, to_bus)
        bus_tree[to_bus] = add_bus_to_tree(bus_tree, to_bus, from_bus)
    end
    return bus_tree
end

function dfs(distance::Int, target_bus_id::Int, bus_tree::Array{Union{Nothing, Array{Int}}};
             neighbour_dict::Dict{Int, Array{Int}}=Dict{Int, Array{Int}}(), parent_bus_id::Int=-1, depth::Int=0)
    if distance <= -1
        return
    end
    
    neighbour_dict[depth] = add_bus_as_neighbour(neighbour_dict, depth, target_bus_id)
    for child_bus_id in bus_tree[target_bus_id]
        if child_bus_id != parent_bus_id
            depth += 1
            dfs(distance - 1, child_bus_id, bus_tree, 
                neighbour_dict=neighbour_dict, 
                parent_bus_id=target_bus_id, depth=depth)
            depth -= 1
        end
    end
    return neighbour_dict
end

function remove_cycled_ids(neighbour_dict::Dict{Int, Array{Int}})
    tmp_arr = []
    for distance in sort(collect(keys(neighbour_dict)))
        for child_bus in neighbour_dict[distance]
            if child_bus ∉ tmp_arr
                tmp_arr = cat(tmp_arr, child_bus, dims=1)
            else
                neighbour_dict[distance] = filter(x -> x ≠ child_bus, neighbour_dict[distance])
            end
        end
    end
    return neighbour_dict
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

function get_children_of_failed_line(casedata::CaseData, failed_line_id::Int;
                                     distance::Int=1, recursive::Bool=true)
    failed_line = get_line_bus_pairs(casedata, failed_line_id, reverse_dict=false)
    bus_tree = get_bus_tree(casedata)

    children_dict = Dict{Tuple{Int, Int, Int}, Array{Tuple{Int,Int}}}()
    for failed_bus in collect(Iterators.flatten(values(failed_line)))
        neighbour_dict = remove_cycled_ids(dfs(distance, failed_bus, bus_tree))
        if recursive
            for depth in setdiff(keys(neighbour_dict), 0)
                children_dict[(failed_line_id, failed_bus, depth)] = [failed_bus < neighbour_dict[depth][idx] ? (failed_bus, neighbour_dict[depth][idx]) : (neighbour_dict[depth][idx], failed_bus) for idx in eachindex(neighbour_dict[depth])]
            end
        else
            children_dict[(failed_line_id, failed_bus, distance)] = [failed_bus < neighbour_dict[distance][idx] ? (failed_bus, neighbour_dict[distance][idx]) : (neighbour_dict[distance][idx], failed_bus) for idx in eachindex(neighbour_dict[distance])]
        end
    end
    return children_dict
end

function remove_already_failed_line(children_dict::Dict{Tuple{Int, Int, Int}, Array{Tuple{Int,Int}}}, failed_line_id::Int)
    failed_line = first(values(get_line_bus_pairs(casedata, failed_line_id, reverse_dict=false)))
    for (key, line_arr) in children_dict
        if failed_line ∈ line_arr
            children_dict[key] = filter(x -> x ≠ failed_line, line_arr)
        end
    end
    return children_dict
end

function match_to_line_id(rev_line_bus_pairs::Dict{Tuple{Int, Int}, Int}, children_dict::Dict{Tuple{Int, Int, Int}, Array{Tuple{Int,Int}}})
    line_id_arr = []
    for lines in collect(Iterators.flatten(values(children_dict)))
        line_id_arr = cat(line_id_arr, rev_line_bus_pairs[lines], dims=1)
    end
    return unique(line_id_arr)
end
    
function get_second_failed_line_id(casedata::CaseData, initial_failed_line_id::Int;
                                   distance::Int=1, recursive::Bool=true)
    children_dict = get_children_of_failed_line(casedata, initial_failed_line_id, distance=distance, recursive=recursive)
    valid_children_dict = remove_already_failed_line(children_dict, initial_failed_line_id)
    rev_line_bus_pairs = get_line_bus_pairs(casedata, reverse_dict=true)
    valid_lines = match_to_line_id(rev_line_bus_pairs, valid_children_dict)
    return valid_lines
end
