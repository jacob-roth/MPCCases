function complete_file_path(file_path::String)
    return !(file_path[end] ∈ Set(['/',"/"])) ? file_path * "/" : file_path
end

# complete the base_files dictionary if not composed
function complete_base_files(base_files::Dict{String, Array}, file_path::String, file_name::String; T::Type=Float64)
    file_exts = (".bus", ".branch", ".gen", ".gencost", ".phys")
    for f_ext in file_exts
      if !haskey(base_files, f_ext)
        base_files[f_ext] = readdlm(complete_file_path(file_path) * file_name * f_ext, T)
      end
    end
    return base_files
  end

function get_y_idx(file_ext::String, P::Bool, Q::Bool)
    if file_ext == ".bus"
        return P & Q ? [3,4] : P ? [3] : Q ? [4] : nothing
    elseif file_ext == ".gen"
        return P & Q ? [2,3] : P ? [2] : Q ? [3] : nothing
    else
        throw(DomainError(file_ext, "file_ext is not properly defined for the number of boolean parameters given."))
    end
end

function get_y_idx(file_ext::String, c2::Bool, c1::Bool, c0::Bool)
    if file_ext == ".gencost"
        return c2 & c1 & c0 ? [5,6,7] : c2 & c1 ? [5,6] : c2 & c0 ? [5,7] : c1 & c0 ? [6,7] : c2 ? [5] : c1 ? [6] : c0 ? [7] : nothing
    else
        throw(DomainError(file_ext, "file_ext is not properly defined for the number of boolean parameters given."))
    end
end

function get_y_idx(file_ext::String, rateA::Bool)
    if file_ext == ".branch" 
        return rateA ? [6] : nothing
    else
        throw(DomainError(file_ext, "file_ext is not properly defined for the number of boolean parameters given."))
    end
end

function get_write_cols_idx(file_ext::String)
    if file_ext ∈ (".bus", ".gen")
        P, Q = (true, true)
        return get_y_idx(file_ext, P, Q)
    elseif file_ext == ".gencost"
        c2, c1, c0 = (true, true, true)
        return get_y_idx(file_ext, c2, c1, c0)
    elseif file_ext == ".branch"
        rateA = true
        return get_y_idx(file_ext, rateA)
    else
        throw(DomainError(file_ext, "file_ext is not properly defined."))
    end
end

function fill_write_file_path(curr_write_file_path::String, read_file_path::String, overwrite_file::Bool, suffix::String)
    filled_write_file_path = complete_file_path(mkpath(
        overwrite_file                  ?   read_file_path  :
        !isempty(curr_write_file_path)  ?   curr_write_file_path :
                                            read_file_path * suffix
                                            ))
    return filled_write_file_path
end

function fill_write_file_name(file_name::String, write_file_name::String, overwrite_file::Bool)
    filled_write_file_name =    overwrite_file              ?   file_name : 
                                !isempty(write_file_name)   ?   write_file_name : 
                                                                file_name
    return filled_write_file_name
end

function match_length(target::Union{Nothing, String, Real, Array{<:Real}}, N::Int)
    return isa(target, Union{Nothing, String, Real}) ? Tuple(repeat([target], N)) : fill(target, N)
end

function cp_remaining_files(src_path::String, dst_path::String, file_name::String)
    file_exts = (".bus", ".gen", ".gencost", ".branch", ".phys")
    src_file_path = complete_file_path(src_path) * file_name
    dst_file_path = complete_file_path(dst_path) * file_name
    for ext in file_exts
        if !isfile(dst_file_path * ext)
            cp(src_file_path * ext, dst_file_path * ext, force=false)
        end
    end
end

function cp_remaining_files(src_path::Union{String, NTuple{N, String}}, dst_path::Union{String, NTuple{N, String}}, file_name::Union{String, NTuple{N, String}}) where {N}
    src_path = isa(src_path, Union{String, Tuple{String}}) ? match_length(convert(String, src_path), N) : src_path
    dst_path = isa(src_path, Union{String, Tuple{String}}) ? match_length(convert(String, dst_path), N) : dst_path
    file_name = isa(file_name, Union{String, Tuple{String}}) ? match_length(convert(String, file_name), N) : file_name
    for idx in 1:N
        cp_remaining_files(src_path[idx], dst_path[idx], file_name[idx])
    end
end

function get_path_dict(cascades_root::String, case_name::String)
    cascades_root = complete_file_path(cascades_root)
    path_dict = Dict{Symbol, String}()
    path_dict[:casedata_path] = cascades_root * "casedata/" * case_name * "/"
    path_dict[:operatingdata_path] = cascades_root * "operatingdata/" * case_name * "/"
    path_dict[:kmcdata_path] = cascades_root * "kmcdata/"
    return path_dict
end

function get_case_arr(cascades_root::String, case_name::String, file_name::String, file_ext::String)
    path_dict = get_path_dict(cascades_root, case_name)
    case_arr = readdlm(path_dict[:casedata_path] * file_name * file_ext)
    return case_arr
end

function add_gaussian_noise(vals::VecOrMat{<:Real}, mean::Real, sd::Real, seed::Union{Nothing, Int})
    rng = isnothing(seed) ? MersenneTwister() : MersenneTwister(seed)
    dims = size(vals)
    gaussian_noise = randn(rng, dims)
    scaled_gaussian_noise = (sd .* gaussian_noise) .+ mean
    return vals + scaled_gaussian_noise
end

# For Pd, c2, c1, c0, if adjustment in adj_arr is negative, use the values from vals instead. Also used in Pg from xbar
function undo_neg_vals(adj_arr::VecOrMat{<:Real}, start_x_idx::Int, y_idx::Union{Nothing, Vector{Int}}, vals::VecOrMat{<:Real})
    if isnothing(y_idx)
        return adj_arr
    else
        vals_length = size(vals, 1)
        subset_adj_arr = adj_arr[start_x_idx : (start_x_idx + vals_length - 1), y_idx]
        discard_neg_vals = subset_adj_arr .* (subset_adj_arr .>= 0) + vals .* (subset_adj_arr .< 0)
        adj_arr[start_x_idx : (start_x_idx + vals_length - 1), y_idx] = discard_neg_vals
        return adj_arr
    end
end

function get_phys(buses::AbstractArray; Dv::T, Mg::T, Dl::T, Dg::T) where T <: AbstractFloat
    phys = MPCCases.Phys[]
    for i in eachindex(buses)
        if buses.bustype[i] == 1
            push!(phys, MPCCases.Phys(buses.bus_i[i], 0.0, Dl, Dv))
        elseif buses.bustype[i] == 2 || buses.bustype[i] == 3
            push!(phys, MPCCases.Phys(buses.bus_i[i], Mg, Dg, 0.0))
        end
    end
    return phys
end

function get_Bs_psd_adjustments(opfdata::OPFData, options::Dict, verb::Bool=false)
    """
    make `B` matrix diagonally dominant
    """
    buses = opfdata.buses
    Bs_adj = zeros(length(buses))
    Y = computeAdmittanceMatrix(opfdata, options)
    B = -imag.(Y)
    for i in 1:size(B,1)
        r = B[i,:]
        rdiag = r[i]
        roffd = sum(r[[(j .!= i) for j in 1:size(B,1)]])
        if abs(roffd) > rdiag
            adj = -(abs(roffd + rdiag) + 1e-7)
            adj_MVA = adj * opfdata.baseMVA
            Bs_adj[i] = adj_MVA
            if verb; println("adjusting for row $i: ", adj_MVA); end
        end
    end
    return Bs_adj
end

function update_case!(casedata::CaseData; buses=nothing, lines=nothing, generators=nothing)
    if !isnothing(buses)
      # build a dictionary between buses ids and their indexes
      casedata.opf.BusIdx .= mapBusIdToIdx(buses)
      casedata.opf.bus_ref .= findall(buses.bustype .== 3)[1]
    end
  
    if !isnothing(lines)
      # set up the FromLines and ToLines for each bus
      buses = casedata.opf.buses
      fl, tl = mapLinesToBuses(buses, lines, casedata.opf.BusIdx)
      casedata.opf.FromLines = fl
      casedata.opf.ToLines = tl
    end
  
    if !isnothing(generators)
      # generators at each bus
      buses = casedata.opf.buses
      casedata.opf.BusGeners = mapGenersToBuses(buses, generators, casedata.opf.BusIdx)
    end
  end
  