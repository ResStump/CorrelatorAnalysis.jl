module CorrelatorAnalysis

import ADerrors as AD
import BDIO
import Dates
import ForwardDiff as FD
import FiniteDifferences as FDiff
import LinearAlgebra as LA
import LsqFit
import PhysicalConstants.CODATA2018 as PCC18
import CairoMakie as CM
using Makie
import Random
import SpecialFunctions as SF
import Statistics as Stats
import Measurements as M
using LaTeXStrings
import PrecompileTools

export
# Types
    FitResult,

# Functions
    read_uwreal,
    write_uwreal,
    export_samples,
    uwone,
    add_mcid_to_parms!,
    err!,
    derivedobs_fd,
    cov,
    posdef_cov,
    GEVP,
    overlaps_Z,
    effective_energy,
    fit_error,
    fit,
    fit_plateau,
    plot_correlator!,
    plot_correlator,
    plot_autocorrelation!,
    plot_autocorrelation,
    plot_effective_energy!,
    plot_effective_energy,
    plot_error_rectangle!,
    plot_error_rectangle,
    plot_herrorline!,
    plot_herrorline,
    plot_model!,
    plot_model,
    plot_overlaps!,
    tomeas,
    uwreal,
    uwreal_array,
    fold_correlator,
    markov_chain,
    bootstrap_to_uwreal,

# Constants
    ħc


include("parms.jl")
include("IO.jl")
include("utils.jl")
include("computations.jl")
include("plot_functions.jl")

# Precompile module
include("precompilation.jl")

end # module CorrelatorAnalysis
