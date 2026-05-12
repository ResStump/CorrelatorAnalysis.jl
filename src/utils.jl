# Overload base functions
Base.abs(a::AD.uwreal) = a*sign(AD.value(a))
Base.length(a::AD.uwreal) = 1
Base.iterate(a::AD.uwreal) = (a, nothing)
Base.iterate(a::AD.uwreal, state) = nothing

for op in (:*, :/)
    @eval function Base.$op(a::AbstractArray{AD.uwreal}, b::AD.uwreal)
        r = similar(a)
        @inbounds for i in eachindex(a)
            r[i] = Base.$op(a[i], b)
        end
        
        return r
    end

    @eval function Base.$op(a::AD.uwreal, b::AbstractArray{AD.uwreal})
        r = similar(b)
        @inbounds for i in eachindex(b)
            r[i] = Base.$op(a, b[i])
        end
        
        return r
    end
end

Base.isapprox(x::AD.uwreal, y::AD.uwreal; kargs...) = 
    isapprox(AD.value(err!(x)), AD.value(err!(y)); kargs...) &&
    isapprox(AD.err(x), AD.err(y); kargs...)
Base.isapprox(x::AbstractArray{AD.uwreal}, y::AbstractArray{AD.uwreal}; kargs...) = 
    isapprox(AD.value.(err!.(x)), AD.value.(err!.(y)); kargs...) &&
    isapprox(AD.err.(x), AD.err.(y); kargs...)

"""
    tomeas(a::AD.uwreal) -> M.measurement

Convert an `AD.uwreal` object to a `Measurements.measurement` object.
"""
tomeas(a::AD.uwreal) = M.measurement(a.mean, err!(a).err)

"""
    uwreal(data::AbstractVector, mcid::String, window=:auto; S=2.0, calc_err=true) -> a::AD.uwreal

Create an `AD.uwreal` object from the input `data` vector. Specify an unique label `mcid`
for the ensemble. The `window` parameter specifies the summation window for the Γ-method
(`window=1` means autocorrelation is neglected). If it's set to `:auto` Ulli Wolff's
automatic windowing procedure with the given parameter `S` (default is 2.0) is used.

See also: `uwreal_array`.
"""
function uwreal(data::AbstractVector, mcid::String, window=:auto; S=2.0, calc_err=true)
    # Add mcid to parms
    add_mcid_to_parms!(mcid, window, S=S)

    a = AD.uwreal(collect(data), mcid)

    if calc_err
        err!(a)
    end

    return a
end

"""
    uwreal_array(data::AbstractArray, mcid::String, window=:auto, mc_dim=:last; S=2.0, calc_err=true) -> uwdata::Array{AD.uwreal}

Create an array of `AD.uwreal` objects from the input `data` array which is assumed to have
the Monte Carlo (MC) time in the last dimension (default) or in dimension `mc_dim`. Specify
an unique label `mcid` for the ensemble. The `window` parameter specifies the summation
window for the Γ-method (`window=1` means autocorrelation is neglected). If it's set to
`:auto` Ulli Wolff's automatic windowing procedure with the given parameter `S`
(default is 2.0) is used.

See also: `uwreal`.
"""
function uwreal_array(data::AbstractArray, mcid::String, window=:auto, mc_dim=:last;
                      S=2.0, calc_err=true)
    # Add mcid to parms
    add_mcid_to_parms!(mcid, window, S=S)

    # Shape of data
    if mc_dim == :last
        mc_dim = ndims(data)
    elseif mc_dim == :first
        mc_dim = 1
    elseif isinteger(mc_dim) && 1 <= mc_dim <= ndims(data)
        mc_dim = Int(mc_dim)
    else
        throw(ArgumentError("mc_dim it not valid."))
    end

    # Allocate array
    uwdata_dims = [dim for dim in 1:ndims(data) if dim!=mc_dim]
    uwdata = Array{AD.uwreal}(undef, size(data)[uwdata_dims])

    # Loop over each dim exept the mc_dim and create uwreal
    for (idx, d) in enumerate(eachslice(data, dims=Tuple(uwdata_dims)))
        uwdata[idx] = AD.uwreal(collect(d), mcid)
    end

    if calc_err
        err!(uwdata)
    end

    return uwdata
end

"""
    fold_correlator(Cₜ::AbstractVector{AD.uwreal}) -> Cₜ_folded::Vector{AD.uwreal}

Fold the correlator `Cₜ` by averaging the entries `Cₜ[i]` and `Cₜ[Nₜ-i]` for `i in 2:Nₜ/2`
where `Nₜ` is the length of `Cₜ`. The entries `Cₜ[1]` and `Cₜ[Nₜ/2+1]` are left unchanged.
"""
function fold_correlator(Cₜ::AbstractVector{AD.uwreal})
    Nₜ = length(Cₜ)
    Cₜ_folded = Cₜ[1:Nₜ÷2+1]
    Cₜ_folded[2:Nₜ÷2] = 0.5*(Cₜ_folded[2:Nₜ÷2] + Cₜ[Nₜ:-1:Nₜ÷2+2])

    return Cₜ_folded
end

"""
    markov_chain(rng, N::Integer, μ::Real, σ::Real, τ::Real) -> uwdata::Vector{AD.uwreal}
    markov_chain(N::Integer, μ::Real, σ::Real, τ::Real) -> uwdata::Vector{AD.uwreal}

Generate a Markov chain of lenght `N` with mean `μ` and error `σ`
(including autocorrelation) and an integrated autocorrelation time `τ`. Optionally provide
a random number generator rng from the `Random` library.
"""
function markov_chain(rng, N::Integer, μ::Real, σ::Real, τ::Real)
    decay = exp(-1/τ)

    # Generate correlated noise with mean 0 and std 1
    noise = Array{Float64}(undef, N)
    noise[1] = randn(rng)
    for i in 2:N
        noise[i] = randn(rng) + decay*noise[i-1]
    end

    # Normalization
    noise *= √(1 - exp(-2/τ))

    # Set error such that error of mean is σ
    noise *= σ * √(N/2τ)

    return μ .+ noise
end
markov_chain(N::Integer, μ::Real, σ::Real, τ::Real) =
    markov_chain(Random.MersenneTwister(), N, μ, σ, τ)

"""
    bootstrap_to_uwreal(mean, samples, mcid) -> AD.uwreal

Convert the bootstrap samples `samples` with mean `mean` to an `AD.uwreal` object with the
same mean and error. Specify an unique label `mcid` for the ensemble.
"""
function bootstrap_to_uwreal(mean, samples, mcid)
    # Scale up error and correct mean
    samples = (samples .- Stats.mean(samples))*√length(samples) .+ mean
    return uwreal(samples, mcid, 1)
end

"""
    add_systematic_error(a_pref, a_alt, α=0.5; label_corr="", labels_uncorr=nothing) -> Vector{AD.uwreal}

Add a systematic uncertainty to a preferred set of observables `a_pref`, estimated from the
difference to a alternative results `a_alt` (all are vectors of `AD.uwreal`). The
correlation parameter `α` sets how correlated the systematic uncertainties are; it must
satisfy `0 < α < 1` (default is 0.5).

The optional string `label_corr` and array of strings `labels_uncorr` are used to label the correlated and uncorrelated systematic ensemble tags respectively.

# Background

Given a preferred result `a_pref[i]` and an alternative result `a_alt[i]` for each data
point, the difference

    Δx[i] = value(a_pref[i]) - value(a_alt[i])

is used as an estimate of the systematic uncertainty at each point. The systematic
covariance matrix is parametrised as:

    cov_sys = α · diag(Δx)² + (1 - α) · Δx Δxᵀ

where the parameter `α ∈ (0, 1)` controls the correlation structure:
- `α → 1`: fully **uncorrelated** systematics (each point shifts independently)
- `α → 0`: fully **correlated** systematics (all points shift together)

# Implementation via ADerrors ensemble tags

Rather than working with the covariance matrix explicitly, the systematic is injected
directly as new `AD.uwreal` sources, so that all subsequent operations (fits, derived
quantities, etc.) propagate the systematic automatically via automatic differentiation.

The decomposition used is:

    a_pref[i] + √(1-α) · Δx[i] · η_shared + √α · Δx[i] · η_i

where `η_shared` is a single shared noise variable (one synthetic "configuration"),
and `η_i` are independent noise variables per data point. This exactly reproduces
`cov_sys` above:

    Cov(i, j) = (1-α)·Δx[i]·Δx[j]   (from shared tag, i ≠ j)
    Var(i)    = (1-α)·Δx[i]² + α·Δx[i]²  = Δx[i]²  (diagonal, as expected)

# Notes
- For a fully correlated systematic (`α → 0`), prefer a small but nonzero value
  such as `α = 1e-6` rather than exactly `0` to avoid zero-variance uncorrelated
  sources.
"""
function add_systematic_error(a_pref::AbstractVector{AD.uwreal},
                              a_alt::AbstractVector{AD.uwreal}, α=0.5;
                              label_corr="", labels_uncorr=nothing)
    n = length(a_pref)
    if length(a_alt) != n
        throw(ArgumentError("a_pref and a_alt must have the same length, got $(n) and "*
                            "$(length(a_alt))."))
    end
    if !(0.0 < α < 1.0)
        throw(ArgumentError("α must satisfy 0 < α < 1, got α = $α."))
    end    
    if isnothing(labels_uncorr)
        labels_uncorr = ["_$i" for i in 1:n]
    else
        if length(labels_uncorr) != n
            throw(ArgumentError("labels_uncorr must have the same length as a_pref, got " *
                                "$(length(labels_uncorr)) and $(n)."))
        end
        labels_uncorr = "_" .* labels_uncorr
    end

    # Systematic magnitude: use only the mean of the difference
    Δa = AD.value.(a_pref .- a_alt)

    # Correlated part: all points share one ensemble tag
    sys_corr = [AD.uwreal([0.0, sqrt(1.0 - α) * Δa[i]], "sys_corr$(label_corr)")
                for i in 1:n]

    # Uncorrelated part: each point gets its own independent ensemble tag
    tags_uncorr = ["sys_uncorr$(label_corr)_$(labels_uncorr[i])" for i in 1:n]
    sys_uncorr = [AD.uwreal([0.0, sqrt(α) * Δa[i]], tags_uncorr[i]) for i in 1:n]

    return a_pref .+ sys_corr .+ sys_uncorr
end