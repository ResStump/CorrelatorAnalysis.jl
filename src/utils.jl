# Overload base functions
Base.abs(a::AD.uwreal) = a*sign(AD.value(a))

for op in (:+, :-, :*, :/, :^, :atan, :hypot)
    @eval function Base.$op(a::Array{AD.uwreal}, b::AD.uwreal)
        r = similar(a)
        @inbounds for i in 1:length(a)
            r[i] = Base.$op(a[i], b)
        end
        
        return r
    end

    @eval function Base.$op(a::AD.uwreal, b::Array{AD.uwreal})
        r = similar(b)
        @inbounds for i in 1:length(b)
            r[i] = Base.$op(a, b[i])
        end
        
        return r
    end
end

function uwreal_array(data::AbstractArray, mcid, window=:auto, mc_dim=:last;
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
        err!.(uwdata)
    end

    return uwdata
end

function fold_correlator(Cₜ::Vector{AD.uwreal})
    Nₜ = length(Cₜ)
    Cₜ_folded = Cₜ[1:Nₜ÷2+1]
    Cₜ_folded[2:Nₜ÷2] = 0.5*(Cₜ_folded[2:Nₜ÷2] + Cₜ[Nₜ:-1:Nₜ÷2+2])

    return Cₜ_folded
end

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