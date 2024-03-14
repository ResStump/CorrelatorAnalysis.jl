mutable struct Parms
    wpm::Dict{String, Vector{Float64}}

    Parms() = new()
end

# Instance of parms
parms = Parms()


function add_mcid_to_parms!(mcid::String, window=:auto; S=2.0)
    # Check if wpm dict is already defined
    if !isdefined(parms, :wpm)
        parms.wpm = Dict{String, Vector{Float64}}()
    end

    if window === :auto
        parms.wpm[mcid] = Float64[-1, S, -1, -1]
    elseif window >= 1
        parms.wpm[mcid] = Float64[window, -1, -1, -1]
    else
        throw(ArgumentError("window has to be :auto or >= 1"))
    end
end

# Useful constants
ħc = (PCC18.c_0*PCC18.ħ/PCC18.e*1e9).val # MeV*fm

