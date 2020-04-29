function complete_file_path(file_path::String)
    return !(file_path[end] ∈ Set(['/',"/"])) ? file_path * "/" : file_path
end

function fill_write_file_path(curr_write_file_path::String, read_file_path::String, overwrite_file::Bool, suffix::String)
    filled_write_file_path = complete_file_path(mkpath(
        overwrite_file                  ?   read_file_path  :
        !isempty(curr_write_file_path)  ?   curr_write_file_path :
                                            read_file_path * suffix))
    return filled_write_file_path
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
