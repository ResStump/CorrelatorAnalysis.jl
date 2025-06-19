module CorrelatorAnalysis

import ADerrors as AD
import BDIO
import Dates
import ForwardDiff as FD
import FiniteDifferences as FDiff
import LinearAlgebra as LA
import LsqFit
import PhysicalConstants.CODATA2018 as PCC18
import Plots
import Random
import SpecialFunctions as SF
import Statistics as Stats
import Measurements as M
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
