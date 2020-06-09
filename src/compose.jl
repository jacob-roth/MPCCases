# Compose Case Columns

# 1
function compose_cols(read_file_path::String, base_file_name::String, aux_file_path::String, file_ext::String; T::Type=Float64)
    read_file_path = complete_file_path(read_file_path)
    base_file = readdlm(read_file_path * base_file_name * file_ext, T)
    aux_cols = readdlm(read_file_path * aux_file_path * file_ext, T)
    write_cols_idx = get_write_cols_idx(file_ext)
    base_file[:, write_cols_idx] = aux_cols
    return base_file
end
