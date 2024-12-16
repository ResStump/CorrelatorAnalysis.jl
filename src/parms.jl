mutable struct Parms
    wpm::Dict{String, Vector{Float64}}
end

# Instance of parms
parms = Parms(Dict{String, Vector{Float64}}())


"""
    add_mcid_to_parms!(mcid::String, window=:auto; S=2.0)

For the ensemble with label `mcid` add a "window parameter" entry to the global window
parameters `parms.wpm`. The `window` argument specifies the summation window for the
Γ-method that is used to compute the error for this ensemble (`window=1` means
autocorrelation is neglected) . If it's set to `:auto` Ulli Wolff's automatic windowing
procedure with the given parameter `S` (default is 2.0) is used.

If the id `mcid` exists already it is overwritten.
"""
function add_mcid_to_parms!(mcid::String, window=:auto; S=2.0)
    if window === :auto
        parms.wpm[mcid] = Float64[-1, S, -1, -1]
    elseif window >= 1
        parms.wpm[mcid] = Float64[window, -1, -1, -1]
    else
        throw(ArgumentError("window has to be :auto or >= 1"))
    end
end

# Useful constants
"""Reduced Plank constant times speed of light in MeV*fm"""
ħc = (PCC18.c_0*PCC18.ħ/PCC18.e*1e9).val # MeV*fm

