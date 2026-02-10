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
    add_mcid_to_parms!,
    bootstrap_to_uwreal,
    cov,
    derivedobs_fd,
    effective_energy,
    err!,
    export_samples,
    fit_error,
    fit_plateau,
    fit,
    fold_correlator,
    GEVP,
    markov_chain,
    overlaps_Z,
    plot_autocorrelation,
    plot_autocorrelation!,
    plot_correlator,
    plot_correlator!,
    plot_effective_energy,
    plot_effective_energy!,
    plot_error_rectangle,
    plot_error_rectangle!,
    plot_herrorline,
    plot_herrorline!,
    plot_model,
    plot_model!,
    plot_overlaps!,
    posdef_cov,
    read_uwreal,
    tomeas,
    uwone,
    uwreal_array,
    uwreal,
    write_uwreal,

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
