function complete_file_path(file_path::String)
    return !(file_path[end] ∈ Set(['/',"/"])) ? file_path * "/" : file_path
end

# complete the base_files dictionary if not composed
function complete_base_files(base_files::Dict{String, Array}, file_path::String, file_name::String)
    file_exts = (".bus", ".branch", ".gen", ".gencost", ".phys")
    for f_ext in file_exts
      if !haskey(base_files, f_ext)
        base_files[f_ext] = readdlm(complete_file_path(file_path) * file_name * f_ext)
      end
    end
    return base_files
end

function fill_write_file_path(write_file_path::String, read_file_path::String, suffix::String)
    filled_write_file_path = complete_file_path(mkpath(
        !isempty(write_file_path) ? write_file_path : read_file_path * suffix))
    return filled_write_file_path
end

function fill_write_file_name(file_name::String, write_file_name::String)
    filled_write_file_name = !isempty(write_file_name) ? write_file_name : file_name
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
  