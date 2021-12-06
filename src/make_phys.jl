function make_phys(read_file_path::String, file_name::String, gen_MDDv::Array{<:Real}, nongen_MDDv::Array{<:Real};
                   overwrite_file::Bool=false, write_file_path::String="")
    bus_file_path = complete_file_path(read_file_path) * file_name * ".bus"
    bus_arr = readdlm(bus_file_path)
    bus_id = bus_arr[:,1]
    bus_type = bus_arr[:,2]
    num_buses = size(bus_arr, 1)
    phys_arr = [bus_type[idx] in Set([2,3]) ? [bus_id[idx], gen_MDDv...] : [bus_id[idx], nongen_MDDv...] for idx in 1:num_buses]
    
    filled_write_file_path = fill_write_file_path(write_file_path, read_file_path, overwrite_file, "/mapped")
    open(filled_write_file_path * file_name * ".phys", "w") do io
        writedlm(io, phys_arr)
    end
end

# gen_MDDv = [0.0531, 0.05, 0.01]
# nongen_MDDv = [0.0, 0.005, 0.01]