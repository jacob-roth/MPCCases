struct Bus
  bus_i::Int
  bustype::Int
  Pd::Float64
  Qd::Float64
  Gs::Float64
  Bs::Float64
  area::Int
  Vm::Float64
  Va::Float64
  baseKV::Float64
  zone::Int
  Vmax::Float64
  Vmin::Float64
end

struct Line
    id::Int
    from::Int
    to::Int
    r::Float64
    x::Float64
    b::Float64
    rateA::Float64
    rateB::Float64
    rateC::Float64
    ratio::Float64 #TAP
    angle::Float64 #SHIFT
    status::Int
    angmin::Float64
    angmax::Float64
end
Line() = Line(0,0,0,0.,0.,0.,0.,0.,0.,0.,0.,0,0.,0.)

mutable struct Gener
  # .gen fields
  bus::Int
  Pg::Float64
  Qg::Float64
  Qmax::Float64
  Qmin::Float64
  Vg::Float64
  mBase::Float64
  status::Int
  Pmax::Float64
  Pmin::Float64
  Pc1::Float64
  Pc2::Float64
  Qc1min::Float64
  Qc1max::Float64
  Qc2min::Float64
  Qc2max::Float64
  ramp_agc::Float64
  # .gencost fields
  gentype::Int
  startup::Float64
  shutdown::Float64
  n::Int64
  coeff::Array
end

mutable struct OPFData
  buses::StructArray{Bus,1}
  lines::StructArray{Line,1}
  generators::StructArray{Gener,1}
  bus_ref::Int
  baseMVA::Float64
  BusIdx::Dict{Int,Int}    #map from bus ID to bus index
  FromLines::Array         #From lines for each bus (Array of Array)
  ToLines::Array           #To lines for each bus (Array of Array)
  BusGenerators::Array     #list of generators for each bus (Array of Array)
end

mutable struct Phys
  ID::Int64
  M::Float64
  D::Float64
  Dv::Float64
end

mutable struct CaseData
  opf::OPFData
  phys::Array{Phys,1}
end

function load_case(completed_base_files::Dict, lineOff::Line=Line(); other::Bool=true, baseMVA::Int=100, T::Type=Float64)
  buses, bus_ref = load_buses(completed_base_files)
  lines = load_branches(completed_base_files, lineOff=lineOff)
  generators = load_generators(completed_base_files, baseMVA)
  phys = load_phys(completed_base_files, other=other)

  # build a dictionary between buses ids and their indexes
  busIdx = mapBusIdToIdx(buses)

  # set up the FromLines and ToLines for each bus
  FromLines,ToLines = mapLinesToBuses(buses, lines, busIdx)

  # generators at each bus
  BusGeners = mapGenersToBuses(buses, generators, busIdx)

  # return opf
  opf = OPFData(StructArray(buses), StructArray(lines), StructArray(generators), bus_ref, baseMVA, busIdx, FromLines, ToLines, BusGeners)
  if other == true
    CD = CaseData(opf, phys)
    return CD
  else
    return opf
  end
end

# For loading cases with composing composite files
function load_case(read_file_path::String, base_file_name::String, aux_file_name::Union{String, NTuple{N, String}}, file_ext::Union{String, NTuple{N, String}}, lineOff::Line=Line(); other::Bool=true, baseMVA::Int=100, T::Type=Float64) where {N}
  base_files = compose_file(read_file_path, base_file_name, aux_file_name, file_ext, T=T)
  completed_base_files = complete_base_files(base_files, read_file_path, base_file_name, T=T)
  return load_case(completed_base_files, lineOff, other=other, baseMVA=baseMVA, T=T)
end

# For loading cases without composing composite files
function load_case(case_name::String, case_path::String, lineOff::Line=Line(); other::Bool=true, baseMVA::Int=100, T::Type=Float64)
  base_files = Dict{String, Array}()
  completed_base_files = complete_base_files(base_files, case_path, case_name, T=T)
  return load_case(completed_base_files, lineOff, other=other, baseMVA=baseMVA, T=T)
end

function load_buses(completed_base_files::Dict)
  bus_arr = completed_base_files[".bus"]
  num_buses = size(bus_arr,1)
  buses = Array{Bus,1}(undef, num_buses)
  bus_ref=-1

  for i in 1:num_buses
    @assert bus_arr[i,1]>0  # don't support nonpositive bus ids
    bus_arr[i,9] *= pi/180  # ANIRUDH: Bus is an immutable struct. Modify bus_arr itself
    buses[i] = Bus(bus_arr[i,1:13]...)
    # buses[i].Va *= pi/180 # ANIRUDH: See previous comment
    if buses[i].bustype==3
      if bus_ref>0
        error("More than one reference bus present in the data")
      else
         bus_ref=i
      end
    end
  end

  return buses, bus_ref
end

function load_branches(completed_base_files::Dict; lineOff::Line=Line())
  branch_arr = completed_base_files[".branch"]
  num_lines = size(branch_arr,1)
  lines_on = findall((branch_arr[:,11].>0) .& ((branch_arr[:,1].!=lineOff.from) .| (branch_arr[:,2].!=lineOff.to)))
  num_on   = length(lines_on)

  if lineOff.from>0 && lineOff.to>0
    println("opf_loaddata: was asked to remove line from,to=", lineOff.from, ",", lineOff.to)
  end

  if length(findall(branch_arr[:,11].==0))>0
    println("opf_loaddata: ", num_lines-length(findall(branch_arr[:,11].>0)), " lines are off and will be discarded (out of ", num_lines, ")")
  end

  lines = Array{Line,1}(undef, num_on)

  lit=0
  for i in lines_on
    @assert branch_arr[i,11] == 1  #should be on since we discarded all other
    lit += 1
    lines[lit] = Line(lit, branch_arr[i, 1:13]...)
    # if lines[lit].angmin>-360 || lines[lit].angmax<360
    #   error("Bounds of voltage angles are still to be implemented.")
    # end
  end
  @assert lit == num_on

  return lines
end

function load_generators(completed_base_files::Dict, baseMVA::Int)
  gen_arr = completed_base_files[".gen"]
  costgen_arr = completed_base_files[".gencost"]
  num_gens = size(gen_arr,1)

  @assert num_gens == size(costgen_arr,1)

  gens_on=findall(gen_arr[:,8] .== 1.0); num_on=length(gens_on)
  if num_gens-num_on>0
    println("loaddata: ", num_gens-num_on, " generators are off and will be discarded (out of ", num_gens, ")")
  end

  generators = Array{Gener,1}(undef, num_on)
  i=0

  for git in gens_on
    i += 1
    generators[i] = Gener(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, Int64[0]) #gen_arr[i,1:end]...)

    generators[i].bus      = gen_arr[git,1]
    generators[i].Pg       = gen_arr[git,2] / baseMVA
    generators[i].Qg       = gen_arr[git,3] / baseMVA
    generators[i].Qmax     = gen_arr[git,4] / baseMVA
    generators[i].Qmin     = gen_arr[git,5] / baseMVA
    generators[i].Vg       = gen_arr[git,6]
    generators[i].mBase    = gen_arr[git,7]
    generators[i].status   = gen_arr[git,8]
    @assert generators[i].status==1
    generators[i].Pmax     = gen_arr[git,9]  / baseMVA
    generators[i].Pmin     = gen_arr[git,10] / baseMVA
    generators[i].Pc1      = gen_arr[git,11]
    generators[i].Pc2      = gen_arr[git,12]
    generators[i].Qc1min   = gen_arr[git,13]
    generators[i].Qc1max   = gen_arr[git,14]
    generators[i].Qc2min   = gen_arr[git,15]
    generators[i].Qc2max   = gen_arr[git,16]
    generators[i].gentype  = costgen_arr[git,1]
    generators[i].startup  = costgen_arr[git,2]
    generators[i].shutdown = costgen_arr[git,3]
    generators[i].n        = costgen_arr[git,4]
    if generators[i].gentype == 1
      generators[i].coeff = costgen_arr[git,5:end]
      error("Piecewise linear costs remains to be implemented.")
    else
      if generators[i].gentype == 2
        generators[i].coeff = costgen_arr[git,5:end]
        #println(generators[i].coeff, " ", length(generators[i].coeff), " ", generators[i].coeff[2])
      else
        generators[i].coeff = costgen_arr[git,5:end]
        # error("Invalid generator cost model in the data.") ## NOTE: JR removed because I need slack/reference...
      end
    end
  end

  return generators
end

function load_phys(completed_base_files::Dict; other::Bool=true)
  if other == true
    phys_arr = completed_base_files[".phys"]
    phys = []
    for i in 1:size(phys_arr,1)
      p = Phys(phys_arr[i,:]...)
      push!(phys, p)
    end
  end
  return phys
end


# Builds a map from lines to buses.
# For each line we store an array with zero or one element containing
# the  'From' and 'To'  bus number.
function mapLinesToBuses(buses, lines, busDict)
  nbus = length(buses)
  FromLines = [Int[] for b in 1:nbus]
  ToLines   = [Int[] for b in 1:nbus]
  for i in 1:length(lines)
    busID = busDict[lines[i].from]
    @assert 1<= busID <= nbus
    push!(FromLines[busID], i)

    busID = busDict[lines[i].to]
    @assert 1<= busID  <= nbus
    push!(ToLines[busID], i)
  end
  return FromLines,ToLines
end

# Builds a mapping between bus ids and bus indexes
#
# Returns a dictionary with bus ids as keys and bus indexes as values
function mapBusIdToIdx(buses)
  dict = Dict{Int,Int}()
  for b in 1:length(buses)
    @assert !haskey(dict,buses[b].bus_i)
    dict[buses[b].bus_i] = b
  end
  return dict
end
function mapIdxToBusId(opfdata::OPFData)
  return Dict(b=>a for (a,b) in opfdata.BusIdx)
end

# Builds a map between buses and generators.
# For each bus we keep an array of corresponding generators number (as array).
#
# (Can be more than one generator per bus)
function mapGenersToBuses(buses, generators,busDict)
  gen2bus = [Int[] for b in 1:length(buses)]
  for g in 1:length(generators)
    busID = busDict[ generators[g].bus ]
    #@assert(0==length(gen2bus[busID])) #at most one generator per bus
    push!(gen2bus[busID], g)
  end
  return gen2bus
end

# convert the above into a dict
function mapGenersToBuses(opfdata::OPFData)
  GEN_idx_to_id = Dict((v,k) for (k,v) in enumerate(mapGenersToBuses(opfdata.buses, opfdata.generators, mapIdxToBusId(opfdata))))
  D = Dict()
  for k in keys(GEN_idx_to_id)
    for kk in k
      push!(D, kk => get(GEN_idx_to_id, k, :NA))
    end
  end
  return D
end
