function computeAdmitances(lines, buses, baseMVA;
    lossless::Bool=false, remove_Bshunt::Bool=false, remove_tap::Bool=false,
    verb::Bool=false, loss_scale::AbstractFloat=1.0, nonneg_Bshunt::Bool=false, nobus_Bshunt::Bool=false)

    """ note:
    (1) `remove_Bshunt` refers to both branch and bus shunts
    (2) `remove_tap` removes line ratio and line angle tap adjustments
    (3) `loss_scale` scales the real components of the admittance matrix by `loss_scale` amount
    """

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
        Ys = 1/((lossless ? 0.0 : loss_scale * lines[i].r) + lines[i].x*im)
        #assign nonzero tap ratio
        if !remove_tap
            tap = (lines[i].ratio == 0) ? (1.0) : (lines[i].ratio)
        else
            tap = 1.0
        end

        #add phase shifters
        if (!lossless && !remove_tap)
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

        if remove_tap
            if !iszero(lines[i].ratio) && verb
                println("warning: lossless assumption changes ratio from ", lines[i].ratio, " to 1 for line ", lines[i].from, " -> ", lines[i].to)
            end
        end
        if lossless
            if !iszero(lines[i].r) && verb
                println("warning: lossless assumption changes r from ", lines[i].r, " to 0 for line ", lines[i].from, " -> ", lines[i].to)
            end
            if !iszero(lines[i].angle) && verb
                println("warning: lossless assumption changes angle from ", lines[i].angle, " to 0 for line ", lines[i].from, " -> ", lines[i].to)
            end
        end
        #@printf("[%4d]  tap=%12.9f   %12.9f\n", i, real(tap), imag(tap));
    end

    nbuses = length(buses)
    YshR = zeros(nbuses)
    YshI = zeros(nbuses)
    for i in 1:nbuses
        YshR[i] = (lossless ? 0.0 : (loss_scale * buses[i].Gs / baseMVA))
        Bs = buses[i].Bs

        if i in 1:10
            YshI[i] = (remove_Bshunt ? 0.0 : (Bs / baseMVA))
        else
            YshI[i] = (remove_Bshunt ? 0.0 : (nobus_Bshunt ? 0.0 : Bs / baseMVA)) ## JR: remove bus shunt
        end
        if lossless && !iszero(buses[i].Gs) && verb
            println("warning: lossless assumption changes Gshunt from ", buses[i].Gs, " to 0 for bus ", i)
        end
        if remove_Bshunt && !iszero(buses[i].Bs) && verb
            println("warning: remove-Bshunt assumption changes Bshunt from ", buses[i].Bs, " to 0 for bus ", i)
        end
        if loss_scale != 1.0 && !iszero(buses[i].Gs) && verb
            println("warning: loss-scale assumption changes Gshunt from ", buses[i].Gs, " to $(loss_scale * buses[i].Gs) for bus ", i)
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

function computeAdmitances(opfdata::OPFData;
lossless::Bool=false, remove_Bshunt::Bool=false, remove_tap::Bool=false, verb::Bool=false, loss_scale::AbstractFloat=1.0, nonneg_Bshunt::Bool=false, nobus_Bshunt::Bool=false)

    return computeAdmitances(opfdata.lines, opfdata.buses, opfdata.baseMVA; 
    lossless=lossless, remove_Bshunt=remove_Bshunt, remove_tap=remove_tap, verb=verb, loss_scale=loss_scale, nonneg_Bshunt=nonneg_Bshunt, nobus_Bshunt=nobus_Bshunt)
end

computeAdmittances = computeAdmitances

function computeAdmittanceMatrix(lines, buses, baseMVA, busDict;
          lossless::Bool=true, remove_Bshunt::Bool=true, remove_tap::Bool=true,
          sparse::Bool=true, verb::Bool=false, loss_scale::AbstractFloat=1.0, nonneg_Bshunt::Bool=false, nobus_Bshunt::Bool=false)

    YffR, YffI, YttR, YttI, YftR, YftI, YtfR, YtfI, YshR, YshI = computeAdmitances(lines, buses, baseMVA;
                                            lossless=lossless, remove_Bshunt=remove_Bshunt, remove_tap=remove_tap,
                                            verb=verb, loss_scale=loss_scale, nonneg_Bshunt=nonneg_Bshunt, nobus_Bshunt=nobus_Bshunt)
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
        return im .* Y
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
    nonneg_Bshunt = haskey(options, :nonneg_Bshunt) ? options[:nonneg_Bshunt] : false
    nobus_Bshunt = haskey(options, :nobus_Bshunt) ? options[:nobus_Bshunt] : false
    remove_tap = haskey(options, :remove_tap) ? options[:remove_tap] : false
    verb = haskey(options, :verb) ? options[:verb] : false
    loss_scale = haskey(options, :loss_scale) ? options[:loss_scale] : 1.0

    if lossless && !current_rating
        println("warning: lossless assumption requires `current_rating` instead of `power_rating`\n")
        current_rating = true
    end

    lines = opfdata.lines; buses = opfdata.buses; generators = opfdata.generators; baseMVA = opfdata.baseMVA
    busIdx = opfdata.BusIdx; FromLines = opfdata.FromLines; ToLines = opfdata.ToLines; BusGeners = opfdata.BusGenerators;
    nbus = length(buses); nline = length(lines); ngen  = length(generators)

    return computeAdmittanceMatrix(lines, buses, baseMVA, busIdx;
    lossless=lossless, remove_Bshunt=remove_Bshunt, remove_tap=remove_tap, verb=verb, loss_scale=loss_scale, nonneg_Bshunt=nonneg_Bshunt, nobus_Bshunt=nobus_Bshunt)
end