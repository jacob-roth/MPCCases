function pad_gen(read_file_path::String, file_name::String; 
                 num_zeros::Int=11, overwrite_file::Bool=false, write_file_path::String="")
    gen_file_path = complete_file_path(read_file_path) * file_name * ".gen"
    gen_arr = readdlm(gen_file_path)
    num_gens = size(gen_arr, 1)
    padding = fill(0.0, num_zeros)
    padded_gen_arr = [[gen_arr[idx,:]..., padding...] for idx in 1:num_gens]
    
    filled_write_file_path = fill_write_file_path(write_file_path, read_file_path, overwrite_file, "/mapped")
    open(filled_write_file_path * file_name * ".gen", "w") do io
        writedlm(io, padded_gen_arr)
    end
end
