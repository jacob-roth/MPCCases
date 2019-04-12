const path = pwd() * "/cases/"
const cases = ["case9","case30"]
# const case = "case5"
const tol = 1e-9

import Pkg
Pkg.activate(dirname(dirname(dirname(path))))
Pkg.instantiate()
using Test
using Printf
using MAT
using SparseArrays, LinearAlgebra
include("../src/MPCCases.jl")

## -----------------------------------------------------------------------------
## load
## -----------------------------------------------------------------------------
@testset "Admittance Matrix" begin
for case in cases
    opfdata = MPCCases.load_case(case, path, other=false);

    # shortcuts for compactness
    lines = opfdata.lines; buses = opfdata.buses; generators = opfdata.generators; baseMVA = opfdata.baseMVA
    busIdx = opfdata.BusIdx; FromLines = opfdata.FromLines; ToLines = opfdata.ToLines; BusGeners = opfdata.BusGenerators;
    nbus = length(buses); nline = length(lines); ngen  = length(generators)

    # branch admitances
    Y = MPCCases.computeAdmittanceMatrix(lines, buses, baseMVA, busIdx; lossless=false, remove_Bshunt=false)
    if case == "case9"
        YY = matread("Y9.mat")
    elseif case == "case30"
        YY = matread("Y30.mat")
    end
    @test norm(YY["Y"] - Y) <= tol
end
end # testset