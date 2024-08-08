module CorrelatorAnalysis

import ADerrors as AD
import LinearAlgebra as LA
import LsqFit
import Plots
import Statistics as Stats
import Random
import PhysicalConstants.CODATA2018 as PCC18
using LaTeXStrings
import PrecompileTools

export ħc, add_mcid_to_parms!
export uwreal_array, fold_correlator, markov_chain
export err!, effective_mass, fit, fit_plateau
export plot_correlator!, plot_correlator, plot_autocorrelation
export plot_effective_mass!, plot_effective_mass
export plot_error_rectangle!, plot_error_rectangle, plot_model!, plot_model

include("parms.jl")
include("utils.jl")
include("computations.jl")
include("plot_functions.jl")

# Precompile module
include("precompilation.jl")

end # module CorrelatorAnalysis
