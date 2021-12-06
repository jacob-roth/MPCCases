# Aggregate the ranks from multiple vectors
# TODO: integrate with kmcsim_sensitivity.jl

# For vectors of idx's, get an array of idx's with vectors as rows of the idx_arr
function get_idx_arr(idx...)
    lengths = map(length, idx)
    @assert all(x -> x == first(lengths), lengths) == true
    idx_arr = Array{Int}(undef, length(idx), first(lengths))
    for i in eachindex(idx)
        idx_arr[i, :] = idx[i]
    end
    return idx_arr
end

# For each idx vector in idx_arr, get the array of ranks of each vector based on the column index number
function get_rank_arr(idx...)
    idx_arr = get_idx_arr(idx...)
    rank_arr = similar(idx_arr)
    for i in 1:size(idx_arr, 1)
        for j in 1:size(idx_arr, 2)
            rank_arr[i,j] = findfirst(isequal(j), idx_arr[i, :])
        end
    end
    return rank_arr
end

# Get the aggregate rank array from the array of ranks by taking the average rank across rows
function get_agg_rank_arr(idx...)
    rank_arr = get_rank_arr(idx...)
    agg_rank_arr = sum(rank_arr, dims=1) ./ size(rank_arr, 1)
    return agg_rank_arr
end

# Get the permutation of the aggregate rank array
function get_agg_rank_perm(idx...)
    agg_rank_arr = get_agg_rank_arr(idx...)
    return sortperm(vec(agg_rank_arr))
end
