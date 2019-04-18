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
Line() = Line(0,0,0.,0.,0.,0.,0.,0.,0.,0.,0,0.,0.)

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

function load_case(case_name, case_path, lineOff=Line(); other::Bool=true)
  if ~(case_path[end] ∈ Set(['/',"/"]))
        case_path = case_path * "/"
  end
  case_name = case_path * case_name

  #
  # load buses
  #
  bus_arr = readdlm(case_name * ".bus")
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
    # println("bus ", i, " ", buses[i].Vmin, "      ", buses[i].Vmax)
  end

  #
  # load branches/lines
  #
  branch_arr = readdlm(case_name * ".branch")
  num_lines = size(branch_arr,1)
  lines_on = findall((branch_arr[:,11].>0) .& ((branch_arr[:,1].!=lineOff.from) .| (branch_arr[:,2].!=lineOff.to)) )
  num_on   = length(lines_on)

  if lineOff.from>0 && lineOff.to>0
    println("opf_loaddata: was asked to remove line from,to=", lineOff.from, ",", lineOff.to)
    #println(lines_on, branch_arr[:,1].!=lineOff.from, branch_arr[:,2].!=lineOff.to)
  end
  if length(findall(branch_arr[:,11].==0))>0
    println("opf_loaddata: ", num_lines-length(findall(branch_arr[:,11].>0)), " lines are off and will be discarded (out of ", num_lines, ")")
  end

  lines = Array{Line,1}(undef, num_on)

  lit=0
  for i in lines_on
    @assert branch_arr[i,11] == 1  #should be on since we discarded all other
    lit += 1
    lines[lit] = Line(branch_arr[i, 1:13]...)
    if lines[lit].angmin>-360 || lines[lit].angmax<360
      error("Bounds of voltage angles are still to be implemented.")
    end
  end
  @assert lit == num_on

  #
  # load generators
  #
  gen_arr = readdlm(case_name * ".gen")
  costgen_arr = readdlm(case_name * ".gencost")
  num_gens = size(gen_arr,1)

  baseMVA=100

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

  # build a dictionary between buses ids and their indexes
  busIdx = mapBusIdToIdx(buses)

  # set up the FromLines and ToLines for each bus
  FromLines,ToLines = mapLinesToBuses(buses, lines, busIdx)

  # generators at each bus
  BusGeners = mapGenersToBuses(buses, generators, busIdx)

  #
  # load physical
  #
  if other == true
    phys_arr = readdlm(case_name * ".phys")
    phys = []
    for i in 1:size(phys_arr,1)
      p = Phys(phys_arr[i,:]...)
      push!(phys, p)
    end
  end

  #
  # return
  #
  opf = OPFData(StructArray(buses), StructArray(lines), StructArray(generators), bus_ref, baseMVA, busIdx, FromLines, ToLines, BusGeners)
  if other == true
    CD = CaseData(opf, phys)
    return CD
  else
    return opf
  end
end

function computeAdmitances(lines, buses, baseMVA; lossless::Bool=false, remove_Bshunt::Bool=false)
  """ note: Bshunt refers to both branch and bus shunts... """
  nlines = length(lines)
  YffR=Array{Float64}(undef, nlines)
  YffI=Array{Float64}(undef, nlines)
  YttR=Array{Float64}(undef, nlines)
  YttI=Array{Float64}(undef, nlines)
  YftR=Array{Float64}(undef, nlines)
  YftI=Array{Float64}(undef, nlines)
  YtfR=Array{Float64}(undef, nlines)
  YtfI=Array{Float64}(undef, nlines)

  for i in 1:nlines
    @assert lines[i].status == 1
    Ys = 1/((lossless ? 0.0 : lines[i].r) + lines[i].x*im)
    #assign nonzero tap ratio
    tap = (lines[i].ratio == 0) ? (1.0) : (lines[i].ratio)

    #add phase shifters
    if (!lossless)
      tap *= exp(lines[i].angle * pi/180 * im)
    end

    Ytt = Ys + (remove_Bshunt ? 0.0 : lines[i].b/2*im)  ## JR: remove branch shunt susceptance
    Yff = Ytt / (tap*conj(tap))
    Yft = -Ys / conj(tap)
    Ytf = -Ys / tap

    #split into real and imag parts
    YffR[i] = real(Yff); YffI[i] = imag(Yff)
    YttR[i] = real(Ytt); YttI[i] = imag(Ytt)
    YtfR[i] = real(Ytf); YtfI[i] = imag(Ytf)
    YftR[i] = real(Yft); YftI[i] = imag(Yft)

    if lossless
      if !iszero(lines[i].r)
        println("warning: lossless assumption changes r from ", lines[i].r, " to 0 for line ", lines[i].from, " -> ", lines[i].to)
      end
      if !iszero(lines[i].angle)
        println("warning: lossless assumption changes angle from ", lines[i].angle, " to 0 for line ", lines[i].from, " -> ", lines[i].to)
      end
    end
    #@printf("[%4d]  tap=%12.9f   %12.9f\n", i, real(tap), imag(tap));
  end

  nbuses = length(buses)
  YshR = zeros(nbuses)
  YshI = zeros(nbuses)
  for i in 1:nbuses
    YshR[i] = (lossless ? 0.0 : (buses[i].Gs / baseMVA))
    YshI[i] = (remove_Bshunt ? 0.0 : (buses[i].Bs / baseMVA)) ## JR: remove bus shunt
    if lossless && !iszero(buses[i].Gs)
      println("warning: lossless assumption changes Gshunt from ", buses[i].Gs, " to 0 for bus ", i)
    end
    if remove_Bshunt && !iszero(buses[i].Bs)
      println("warning: remove-Bshunt assumption changes Bshunt from ", buses[i].Bs, " to 0 for bus ", i)
    end
  end

  @assert 0==length(findall(isnan.(YffR)))+length(findall(isinf.(YffR)))
  @assert 0==length(findall(isnan.(YffI)))+length(findall(isinf.(YffI)))
  @assert 0==length(findall(isnan.(YttR)))+length(findall(isinf.(YttR)))
  @assert 0==length(findall(isnan.(YttI)))+length(findall(isinf.(YttI)))
  @assert 0==length(findall(isnan.(YftR)))+length(findall(isinf.(YftR)))
  @assert 0==length(findall(isnan.(YftI)))+length(findall(isinf.(YftI)))
  @assert 0==length(findall(isnan.(YtfR)))+length(findall(isinf.(YtfR)))
  @assert 0==length(findall(isnan.(YtfI)))+length(findall(isinf.(YtfI)))
  @assert 0==length(findall(isnan.(YshR)))+length(findall(isinf.(YshR)))
  @assert 0==length(findall(isnan.(YshI)))+length(findall(isinf.(YshI)))
  if lossless
    @assert 0==length(findall(!iszero, YffR))
    @assert 0==length(findall(!iszero, YttR))
    @assert 0==length(findall(!iszero, YftR))
    @assert 0==length(findall(!iszero, YtfR))
    @assert 0==length(findall(!iszero, YshR))
  end

  return YffR, YffI, YttR, YttI, YftR, YftI, YtfR, YtfI, YshR, YshI
end
computeAdmittances(lines, buses, baseMVA; lossless::Bool=false, remove_Bshunt::Bool=false) = computeAdmitances(lines, buses, baseMVA; lossless::Bool=false, remove_Bshunt::Bool=false)
function computeAdmittances(opfdata::OPFData; lossless::Bool=false, remove_Bshunt::Bool=false)
  return computeAdmittances(opfdata.lines, opfdata.buses, opfdata.baseMVA; lossless=lossless, remove_Bshunt=remove_Bshunt)
end
computeAdmitances(opfdata::OPFData; lossless::Bool=false, remove_Bshunt::Bool=false) = computeAdmittances(opfdata::OPFData; lossless::Bool=false, remove_Bshunt::Bool=false)

function computeAdmittanceMatrix(lines, buses, baseMVA, busDict; lossless::Bool=true, remove_Bshunt::Bool=true, sparse::Bool=true)
  YffR, YffI, YttR, YttI, YftR, YftI, YtfR, YtfI, YshR, YshI = computeAdmitances(lines, buses, baseMVA; lossless=lossless, remove_Bshunt=remove_Bshunt)
  nbuses = length(buses)
  nlines = length(lines)

  if lossless
    Y = zeros(Float64, nbuses, nbuses)
    for l in 1:nlines
      i = busDict[lines[l].from]
      j = busDict[lines[l].to]
      Y[i,j] += YftI[l]
      Y[j,i] += YtfI[l]
      Y[i,i] += YffI[l]
      Y[j,j] += YttI[l]
    end
    if remove_Bshunt == false
      for i in 1:nbuses
        Y[i,i] += YshI[i]
      end
    end
    return Y
  else
    if sparse
      B = spzeros(Float64, nbuses, nbuses)
      G = spzeros(Float64, nbuses, nbuses)
    else
      B = zeros(Float64, nbuses, nbuses)
      G = zeros(Float64, nbuses, nbuses)
    end
    for l in 1:nlines
      i = busDict[lines[l].from]
      j = busDict[lines[l].to]
      B[i,j] += YftI[l]
      B[j,i] += YtfI[l]
      B[i,i] += YffI[l]
      B[j,j] += YttI[l]

      G[i,j] += YftR[l]
      G[j,i] += YtfR[l]
      G[i,i] += YffR[l]
      G[j,j] += YttR[l]
    end
    if remove_Bshunt == false
      for i in 1:nbuses
        B[i,i] += YshI[i]
      end
    end
    for i in 1:nbuses
      G[i,i] += YshR[i]
    end
    return G + im*B
  end
end
function computeAdmittanceMatrix(opfdata::OPFData, options::Dict=Dict())
  # parse options
  lossless = haskey(options, :lossless) ? options[:lossless] : false
  current_rating = haskey(options, :current_rating) ? options[:current_rating] : false
  remove_Bshunt = haskey(options, :remove_Bshunt) ? options[:remove_Bshunt] : false
  if lossless && !current_rating
      println("warning: lossless assumption requires `current_rating` instead of `power_rating`\n")
      current_rating = true
  end
  lines = opfdata.lines; buses = opfdata.buses; generators = opfdata.generators; baseMVA = opfdata.baseMVA
  busIdx = opfdata.BusIdx; FromLines = opfdata.FromLines; ToLines = opfdata.ToLines; BusGeners = opfdata.BusGenerators;
  nbus = length(buses); nline = length(lines); ngen  = length(generators)
  return computeAdmittanceMatrix(lines, buses, baseMVA, busIdx; lossless=lossless, remove_Bshunt=remove_Bshunt)
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
