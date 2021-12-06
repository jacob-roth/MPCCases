function apply_mapping_to_file(read_file_path::String, file_name::String, file_ext::String; 
                               overwrite_file::Bool=false, write_file_path::String="")
    mapping = get_id_mapping(read_file_path, file_name, ".bus")
    full_file_path = complete_file_path(read_file_path) * file_name * file_ext
    arr = readdlm(full_file_path)
    num_rows = size(arr, 1)
    if file_ext in Set([".bus", ".gen", ".phys"])
        mapped_arr = [mapping[arr[idx, 1]] for idx in 1:num_rows]
        arr[:,1] = mapped_arr
    elseif file_ext == ".branch"
        fbus_mapped_arr = [mapping[arr[idx, 1]] for idx in 1:num_rows]
        tbus_mapped_arr = [mapping[arr[idx, 2]] for idx in 1:num_rows]
        arr[:, 1] = fbus_mapped_arr
        arr[:, 2] = tbus_mapped_arr
    else
        error("No bus_id mapping to apply to this file extension")
    end

    filled_write_file_path = fill_write_file_path(write_file_path, read_file_path, overwrite_file, "/mapped")
    open(filled_write_file_path * file_name * file_ext, "w") do io
        writedlm(io, arr)
    end
end

function apply_mapping_to_file(read_file_path::String, file_name::String, file_exts:: Tuple{String, Vararg{String}}; 
                               overwrite_file::Bool=false, write_file_path::String="")
    for file_ext in file_exts
        apply_mapping_to_file(read_file_path, file_name, file_ext, overwrite_file=overwrite_file, write_file_path=write_file_path)
    end
end

# Helper Functions for Apply Mapping to File

function get_bus_ids(read_file_path::String, file_name::String, file_ext::String)
    read_file_path = complete_file_path(read_file_path)
    bus_arr = (file_ext == ".bus") ? readdlm(read_file_path * file_name * file_ext) : error("Not a bus file")
    bus_ids = bus_arr[:, 1]
    @assert allunique(bus_ids)
    return bus_ids
end

function get_id_mapping(read_file_path::String, file_name::String, file_ext::String)
    bus_ids = get_bus_ids(read_file_path, file_name, file_ext)
    num_bus_ids = length(bus_ids)
    ordered_ids = collect(1:num_bus_ids)
    mapping = Dict(bus_ids .=> ordered_ids)
    return mapping
end
