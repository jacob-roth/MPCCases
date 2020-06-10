function complete_file_path(file_path::String)
    return !(file_path[end] ∈ Set(['/',"/"])) ? file_path * "/" : file_path
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

function get_y_idx(file_ext::String, rateA::Bool)
    if (file_ext == ".branch") & rateA
        return 6:6
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

function match_length(target::Union{Nothing, String, Real, Array{Real}}, N::Int)
    return isa(target, Union{Nothing, String, Real}) ? Tuple(repeat([target], N)) : fill(target, N)
end

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

function cp_remaining_files(src_path::Union{String, NTuple{N, String}}, dst_path::Union{String, NTuple{N, String}}, file_name::Union{String, NTuple{N, String}}) where {N}
    src_path = isa(src_path, Union{String, Tuple{String}}) ? match_length(convert(String, src_path), N) : src_path
    dst_path = isa(src_path, Union{String, Tuple{String}}) ? match_length(convert(String, dst_path), N) : dst_path
    file_name = isa(file_name, Union{String, Tuple{String}}) ? match_length(convert(String, file_name), N) : file_name
    for idx in 1:N
        cp_remaining_files(src_path[idx], dst_path[idx], file_name[idx])
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
  