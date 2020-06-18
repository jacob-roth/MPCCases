# Compose Case Columns

# 1: Compose file for single file_ext
function compose_file(read_file_path::String, base_file_name::String, aux_file_name::String, file_ext::String; T::Type=Float64)
    read_file_path = complete_file_path(read_file_path)
    base_file = readdlm(read_file_path * base_file_name * file_ext, T)
    aux_cols = readdlm(read_file_path * aux_file_name * file_ext, T)
    write_cols_idx = get_write_cols_idx(file_ext)
    base_file[:, write_cols_idx] = aux_cols
    return base_file
end

# 2: Compose files for multiple file_exts
function compose_file(read_file_path::String, base_file_name::String, aux_file_name::NTuple{N, String}, file_ext::NTuple{N,String}; T::Type=Float64) where {N}
    @assert file_ext == Tuple(unique(file_ext))
    base_files = Dict{String, Array}()
    for idx in 1:N
        base_file = compose_file(read_file_path, base_file_name, aux_file_name[idx], file_ext[idx])
        base_files[file_ext[idx]] = base_file
    end
    return base_files
end
