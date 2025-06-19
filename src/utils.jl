# Overload base functions
Base.abs(a::AD.uwreal) = a*sign(AD.value(a))

for op in (:+, :-, :*, :/, :^, :atan, :hypot)
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
tomeas(a::AD.uwreal) = M.measurement(a.mean, CA.err!(a).err)

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