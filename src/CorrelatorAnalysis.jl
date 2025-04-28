module CorrelatorAnalysis

import ADerrors as AD
import BDIO
import Dates
import ForwardDiff as FD
import LinearAlgebra as LA
import LsqFit
import PhysicalConstants.CODATA2018 as PCC18
import Plots
import CairoMakie as CM
using Makie
import Random
import SpecialFunctions as SF
import Statistics as Stats
using LaTeXStrings
import PrecompileTools

include("parms.jl")
include("IO.jl")
include("utils.jl")
include("computations.jl")
include("plot_functions.jl")

# Precompile module
include("precompilation.jl")

end # module CorrelatorAnalysis
