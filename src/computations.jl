"""
    err!(a::AD.uwreal) -> a::AD.uwreal

Compute the error of the `AD.uwreal` object `a` in place using the window parameters that
were specified for the ensembles that `a` depends on. For convenience, `a` is also returned.

Look at the documentation of `AD.uwerr` for further informaiton.
"""
function err!(a::AD.uwreal)
    try
        AD.uwerr(a, parms.wpm)
    catch e
        println("Failed to calculate error.")
        println(e)
    end

    return a
end

"""
    err!(a::AbstractArray{AD.uwreal}) -> a::Array{AD.uwreal}

Compute the error of each `AD.uwreal` object in `a` in place using the window parameters
that were specified for the ensembles that `a` depends on. For convenience, `a` is also
returned.

Look at the documentation of `AD.uwerr` for further informaiton.
"""
function err!(a::AbstractArray{AD.uwreal})
    for i in eachindex(a)
        try
            AD.uwerr(a[i], parms.wpm)
        catch e
            println("Failed to calculate error.")
            println(e)
        end
    end

    return a
end

"""
    cov(a::AbstractVector{AD.uwreal}) -> C::Matrix{Float64}

Compute the covariance matrix of the vector `a` using the window parameters that were
specified for the ensembles that `a` depends on.
"""
function cov(a::AbstractVector{AD.uwreal})
    # Copy `a` to avoid changing its error
    a_ = deepcopy(a)

    # Compute covariance matrix using the set window parameters wpm
    cov = AD.cov(a_, parms.wpm)

    return cov
end

"""
    posdef_cov(a::AbstractVector{AD.uwreal}; correlation=false) -> C::Matrix{Float64}

Compute an approximate version of the covariance matrix of the vector `a` which is positive
(semi-)definite (as long as all values in `a` are finite).

If `correlation=true`, compute the correlation matrix instead of the covariance matrix 
(default is `false`).

### Method
This algorithm approximates the integrated autocorrelation time of the offdiagonal covariance
matrix entries `τ_ij` as the geometric mean of those of the observables `a[i]` and
`a[j]` (i.e. `τ_ij = √(τ_i*τ_j)`).

To compute the `τ_i`'s, the largest windows used in `a` to compute
its error (i.e. for each ensemble separately) are used.

This approximation might be bad if `a` depends on multiple ensembles with different
autocorrelation times.

### Acknowledgment
This function is strongly inspired by `pyerrors` (arXiv:hep-lat/0306017) function
obs.covariance.
"""
function posdef_cov(a::AbstractVector{AD.uwreal}; correlation=false)
    # Compute error of `a`
    err!(a)

    # Get biggest windows that are used in `a`
    ids = AD.unique_ids_multi(a, AD.wsg)
    nid = length(ids)
    iw  = zeros(Int64, nid)

    for k in eachindex(a)
        for j in eachindex(a[k].ids)
            for i in 1:nid
                if (a[k].ids[j] == ids[i])
                    iw[i] = max(iw[i], a[k].cfd[j].iw)
                end
            end
        end
    end
    
    wopt = Dict{Int64,Vector{Float64}}()
    for i in 1:nid
        wopt[ids[i]] = Float64[iw[i], -1.0, -1.0, -1.0]
    end

    # Compute the errors in `a` with those windows
    # (Copy `a` to avoid changing its error)
    a_ = deepcopy(a)
    AD.uwerr.(a_, [wopt])
    a_err = AD.err.(a_)

    # Generate window parameters with window 1 (i.e. no autocorrelation)
    wpm0 = deepcopy(wopt)
    for k in keys(wpm0)
        wpm0[k] = [1.0, -1.0, -1.0, -1.0]
    end

    # Get correlation matrix without autocorrelation
    cov0 = AD.cov(a_, wpm0)

    # Scale entries to include autocorrelation
    if correlation
        autocorr_correction = 1 ./ .√LA.diag(cov0)
    else
        autocorr_correction = a_err ./ .√LA.diag(cov0)
    end

    C = LA.diagm(autocorr_correction) * cov0 * LA.diagm(autocorr_correction)
    
    return 0.5*(C + C')
end

function eigvals_AD!(λ_t::AbstractVector{AD.uwreal}, Cₜ::AbstractArray{AD.uwreal, 3},
                     iₜ, i_t₀)
    N_op = size(Cₜ, 2)

    # Compute eigenvectors (sorted in decreasing order of real part of eigenvalue)
    λ, v = LA.eigen(AD.value.(Cₜ[iₜ, :, :]), AD.value.(Cₜ[i_t₀, :, :]),
                    sortby=(λ->-real(λ)))

    # Propagate error to eigenvalues
    for i in 1:N_op
        der = real(conj(v[:, i])*transpose(v[:, i]))
        λ_t[i] = AD.addobs(vec(Cₜ[iₜ, :, :]), vec(der), λ[i])
    end
end

"""
    GEVP(Cₜ::AbstractArray{AD.uwreal, 3}, t₀::Union{Int, Symbol}=:ceil_t_half) -> E_eff::Vector{Vector{{AD.uwreal}}

Solve the Generalized Eigenvalue Problem (GEVP) `Cₜ(t)vₙ = λₙ(t, t₀)Cₜ(t₀)vₙ` and return the
effective energy `Eₙ_eff(t, t₀) = log(λₙ(t, t₀)/λₙ(t+1, t₀))` as a vector of `AD.uwreal`
vectors `E_eff`. The outer index is the eigenvalue index `ₙ` and the inner index is the time
`t`. \\
The parameter `t₀` specifies `t₀` in `Eₙ_eff(t, t₀)`. The options are:
- `:ceil_t_half` which sets `t₀ = ceil(t/2)` (default).
- an `Int`. In that case `t₀` is always the same and the entries for `Eₙ_eff(t, t₀)` with
    `t<t₀` are set to NaN.
"""
function GEVP(Cₜ::AbstractArray{AD.uwreal, 3}, t₀::Union{Int, Symbol}=:ceil_t_half)
    # Get Nₜ and number of operators
    Nₜ, N_op, _ = size(Cₜ)

    if t₀ isa Int
        if t₀ < 0 || t₀ >= Nₜ
            throw(ArgumentError("t₀ must be in the range 0 <= t₀ < Nₜ."))
        end

        # Index of t₀
        i_t₀ = t₀+1

        λ_t = Vector{AD.uwreal}(undef, N_op)
        λ = Array{AD.uwreal, 2}(undef, Nₜ, N_op)
        for iₜ in 1:Nₜ
            if iₜ < i_t₀
                λ[iₜ, :] .= [AD.uwreal(NaN)]
            else
                # Compute eigenvalues
                eigvals_AD!(@view(λ[iₜ, :]), Cₜ, iₜ, i_t₀)
            end
        end

        E_eff = effective_energy.(eachcol(λ), :log)
    elseif t₀ == :ceil_t_half
        λ_t = Vector{AD.uwreal}(undef, N_op)
        λ_tp1 = Vector{AD.uwreal}(undef, N_op)
        E_eff = Array{AD.uwreal, 2}(undef, Nₜ, N_op)
        for iₜ in 1:Nₜ-1
            # Compute t₀ and its index
            t₀ = ceil(Int, iₜ/2)
            i_t₀ = t₀+1

            # Compute eigenvalues
            eigvals_AD!(λ_t, Cₜ, iₜ, i_t₀)
            eigvals_AD!(λ_tp1, Cₜ, iₜ+1, i_t₀)

            # Compute effective energy (set all nonpositive values to NaN)
            ratio = λ_t./λ_tp1

            for i in eachindex(ratio)
                if ratio[i] ≤ 0
                    ratio[i] = AD.uwreal(NaN)
                end
            end
            E_eff[iₜ, :] = @. log(ratio)
        end

        E_eff[end, :] .= [AD.uwreal(NaN)]
        E_eff = eachcol(E_eff)
    else
        throw(ArgumentError("unknown t₀. Use an integer or :ceil_t_half."))
    end
    
    return E_eff
end

"""
    overlaps_Z(Cₜ::Union{AbstractArray{<:Real, 3}, AbstractArray{<:AD.uwreal, 3}}, E_arr::Union{AbstractVector{<:Real}, AbstractVector{<:AD.uwreal}}, t::Int; t₀::Union{Int, Symbol}=:ceil_t_half, normalization::Symbol=:N_inf) -> Z_in::Array{Float64, 2}

Compute the matrix of overlaps `Z_in = <Ω|Oᵢ|n>` of the operator `Oᵢ` with the `n'th`
eigenstates of the Hamiltonian. For that use the correlator matrix
`Cₜ[t+1, i, j] = <Ω|Oᵢ(t)Oⱼ(0)^†|Ω>` and the energies `E_arr` of the lowest eigenstates.
`Z_in` is computed for each operator `Oᵢ` and each eigenstate `n` for which the energy is
given in `E_arr`.

The time `t` and `t₀` is where the generalized eigenvalue problem is solved. The options
for `t₀` are:
- `:ceil_t_half` which sets `t₀ = ceil(t/2)` (default).
- an `Int` with `t₀<t`.
Additionally, specify if and how to normalize the operators `Oᵢ` (and thus `Z_in`).
The options are
- `:N_inf`: normalize `Z_in` such that `Σ_{n=1}^∞ Z_in = Cₜ[1, i, i] = 1` (default).
- `:N_max`: normalize `Z_in` such that `sum(Z_in, dims=2) = 1`
    (equivalent to `Σ_{n=1}^N_max Z_in = 1` for `N_max = length(E_arr)`).
- `:unnormalized`: do not normalize `Z_in`.

# Reference
To compute `Z_in` formula (3.1) in https://doi.org/10.1016/j.nuclphysb.2016.07.024 for
`t_d = t` is used.
"""
function overlaps_Z(Cₜ::Union{AbstractArray{<:Real, 3}, AbstractArray{<:AD.uwreal, 3}},
                   E_arr::Union{AbstractVector{<:Real}, AbstractVector{<:AD.uwreal}},
                   t::Int; t₀::Union{Int, Symbol}=:ceil_t_half,
                   normalization::Symbol=:N_inf)
    # Convert Cₜ and E_arr to Float if necessary
    (Cₜ[1] isa AD.uwreal) && (Cₜ = AD.value.(Cₜ))
    (E_arr[1] isa AD.uwreal) && (E_arr = AD.value.(E_arr))

    # Get time indices
    if t₀ == :ceil_t_half
        t₀ = ceil(Int, t/2)
    elseif t₀ isa Int
        if t₀ >= t
            throw(ArgumentError("t₀ must be smaller than t."))
        end
    else
        throw(ArgumentError("unknown t₀. Use an integer or :ceil_t_half."))
    end
    iₜ, i_t₀ = t+1, t₀+1

    if E_arr[1:end-1] > E_arr[2:end]
        throw(ArgumentError("E_arr must be sorted in ascending order."))
    end

    # Number of levels and operators
    N_levels = length(E_arr)
    N_op = size(Cₜ, 2)

    # Compute generealized eigenvalues and eigenvectors
    # (sorted in decreasing order of real part of eigenvalue)
    λ, v = LA.eigen(Cₜ[iₜ, :, :], Cₜ[i_t₀, :, :], sortby=(λ->-real(λ)))

    Z_in = Array{Float64}(undef, N_op, N_levels)
    for i_n in 1:N_levels
        Z_in[:, i_n] = abs2.(exp(E_arr[i_n]*t/2)/√λ[i_n] * Cₜ[iₜ, :, :]*v[:, i_n])
    end

    if normalization == :N_inf
        return Z_in ./ LA.diag(Cₜ[1, :, :])
    elseif normalization == :N_max
        return stack(eachrow(Z_in) ./ sum.(eachrow(Z_in)), dims=1)
    elseif normalization == :unnormalized
        return Z_in
    else
        throw(ArgumentError("unknown normalization. Use :N_inf, :N_max or :unnormalized."))
    end
end

"""
    effective_energy(Cₜ::AbstractVector{AD.uwreal}, variant=:log; guess=1.0, folded=false) -> E_eff::Vector{AD.uwreal}

Compute the effective energy of the vector `Cₜ` using the specified `variant`.

Optionally, provide the parameter `guess` as an initial value for the root finding
algorithm. If `folded=true`, the correlator is assumed to be folded (only relevant for 
`variant=:cosh` or `variant=:sinh` to determine `Nₜ`).

### Variants
- log: Use the standard effective energy `log(Cₜ(t)/Cₜ(t+1))`.
- cosh: Use the periodicity of the correlator by solving \\
  `Cₜ(t)/Cₜ(t+1) = cosh(E*(t - Nₜ/2)) / cosh(E*(t + 1 - Nₜ/2))` \\
  for E.
- sinh: Use the anti-periodicity of the correlator by solving \\
  `Cₜ(t)/Cₜ(t+1) = sinh(E*(t - Nₜ/2)) / sinh(E*(t + 1 - Nₜ/2))` \\
  for E.
"""
function effective_energy(Cₜ::AbstractVector{AD.uwreal}, variant=:log; guess=1.0,
                          folded=false)
    # Compute error
    err!.(Cₜ)

    if variant in [:exp, :log]
        ratio = Cₜ./circshift(Cₜ, -1)
        
        # Set all nonpositive values to NaN
        for i in eachindex(ratio)
            if ratio[i] ≤ 0
                ratio[i] = AD.uwreal(NaN)
            end
        end

        E_eff = log.(ratio)
    elseif variant in [:cosh, :sinh]
        N_elements = length(Cₜ)
        if folded
            Nₜ = (N_elements - 1) * 2
        else
            Nₜ = N_elements
        end
        
        E_eff = Vector{AD.uwreal}(undef, N_elements)

        # Define fit model
        if variant == :cosh
            f = (E, p) -> cosh(E*(p[1] - Nₜ/2)) / cosh(E*(p[1]+1 - Nₜ/2)) - p[2]
        else
            f = (E, p) -> sinh(E*(p[1] - Nₜ/2)) / sinh(E*(p[1]+1 - Nₜ/2)) - p[2]
        end

        # Loop over all times `t` and compute effective energy
        for (iₜ, t) in enumerate(0:N_elements-1)
            p = [AD.uwreal(Float64(t)), Cₜ[iₜ]/Cₜ[mod1(iₜ+1, N_elements)]]
            # Try to find root, otherwise set value to NaN
            try
                E_eff[iₜ] = abs(AD.root_error(f, Float64(guess), p))
            catch e
                E_eff[iₜ] = AD.uwreal(NaN)
                println("For t = $t (iₜ=$iₜ): ", e)
            end
        end
    else
        throw(ArgumentError("unknown variant to compute effective energy."))
    end

    return E_eff
end

struct FitResult{C, P}
    param::Vector{AD.uwreal}
    χ²_obs::Float64
    χ²_red::C
    dof::Int
    p_value::P
end

"""
    p_value(χ²::Function, data::AbstractVector{AD.uwreal}, p::AbstractArray{<:Real}, W::AbstractMatrix{<:Real}, χ²_obs::Real, dof::Integer; fit_type=nothing, p_value_type=nothing, N_mc=10^5, rng=Random.GLOBAL_RNG) -> p_val

Compute the p-value for the chi-square function `χ²`, the data `data` and the
optimized fit parameters `p`.

### Arguments
- `χ²(p::Vector, d:Vector)`: Function of the parameters `p` and the data `d`. The function is expected to have the following form \\
  `χ²(p, d) = sum_{ij} [d_i - f_i(p)]W_ij[d_j - f_j(p)]` \\
  for the model f_i(p).
- `data::AbstractVector{AD.uwreal}`: The data for which the `χ²` function was
  optimized.
- `p::AbstractArray{<:Real}`: The optimized fit parameters.
- `W::AbstractMatrix{<:Real}`: The weight matrix in `χ²`.
- `χ²_obs::Real`: The χ² value for the optimized parameters.
- `dof::Integer`: The number of degrees of freedom.
- `p_value_type`: Sets how the p-value is computed. It's either `:correlated` or `:general`.
  The first case applies to a correlated fit, then the p-value can be computed using an
  upper incomplete gamma function. For a fit with a different weight matrix `W`
  `p_value_type=:general` should be used. Then the p-value is computed using Monte Carlo
  integration. Default is `nothing`. Either `fit_type` or `p_value_type` have to be given.
- `fit_type`: The type of fit used. Possible choises are `:correlated`, `:correlated_posdef`
  and `:uncorrelated`. The first has the same effect as setting `p_value_type=:correlated`.
  The other two correspond to `p_value_type=:general`. Default is `nothing`.
  Either `fit_type` or `p_value_type` have to be given.
- `N_mc`: The number of Monte Carlo samples used for the integration. Default is `10^5`.
- `rng`: The random number generator. Default is `Random.GLOBAL_RNG`.

### Method
The p-value is computed as described in arXiv:2209.14188 for correlated and uncorrelated
fits. \\
The function `fit` allows to use the approximated covariance matrix computed with 
`posdef_cov`. In that case `W` is treated as a general weight matrix to correct for this.
"""
function p_value(χ²::Function, data::AbstractVector{AD.uwreal}, p::AbstractArray{<:Real},
                 W::AbstractMatrix{<:Real}, χ²_obs::Real, dof::Integer;
                 p_value_type=nothing,  fit_type=nothing, N_mc=10^5, rng=Random.GLOBAL_RNG)
    # Set how to compute p-value, first based on `p_val_type` then on `fit_type`
    if isnothing(p_value_type)
        if fit_type in [:uncorrelated, :correlated_posdef]
            p_value_type = :general
        elseif fit_type == :correlated
            p_value_type = :correlated
        else
            throw(ArgumentError("unknown fit type. Use :uncorrelated, :correlated or "*
                                ":correlated_posdef."))
        end
    elseif !(p_value_type in [:correlated, :general])
        throw(ArgumentError("unknown p-value type. Use :correlated or :general."))
    end

    if p_value_type == :correlated
        # Assumes that the weight matrix `W` is the inverse of the covariance matrix
        p_val = SF.gamma(dof/2, χ²_obs/2)/SF.gamma(dof/2)
    elseif p_value_type == :general
        # Here follow notation in ADerrors

        # Compute Hessian matrix of χ²
        n = length(p)   # Number of fit parameters
        m = length(data) # Number of data

        xav = zeros(Float64, n+m)
        for i in 1:n
            xav[i] = p[i]
        end
        for i in n+1:n+m
            xav[i] = data[i-n].mean
        end
        ccsq(x::Vector) = χ²(view(x, 1:n), view(x, n+1:n+m)) 
        if (n+m < 4)
            cfg = FD.HessianConfig(ccsq, xav, FD.Chunk{1}());
        else
            cfg = FD.HessianConfig(ccsq, xav, FD.Chunk{4}());
        end
            
        hess = Array{Float64}(undef, n+m, n+m)
        FD.hessian!(hess, ccsq, xav, cfg)

        # Compute ν matrix (see arXiv:2209.14188)
        Lm = LA.cholesky(LA.Symmetric(W))
        Li = LA.inv(Lm.L)
        
        hm = view(hess, 1:n, n+1:n+m)
        sm = hm * Li'
        
        maux = sm * sm'
        hi = LA.pinv(maux)
        Pw = hm' * hi * hm

        C = cov(data)
        sqrtC = √(C)
        ν = sqrtC * (W - Pw) * sqrtC

        # Positive eigenvalues of ν
        λ = LA.eigvals(ν)
        if !(λ ≈ real.(λ))
            println("Warning: The eigenvalues of ν are not real and the p-value might "*
                    "therefore be wrong.")
        end
        λ = real.(λ)
        λ_pos = λ[λ .> 0]
        N_ν = length(λ_pos)

        # Compute integral over θ using MC with normaly distributed random numbers
        θ(z) = Float64(sum(@. λ_pos*z^2) ≥ χ²_obs)
        z_arr = randn(rng, (N_ν, N_mc))
        p_val = sum(θ.(eachcol(z_arr)))/N_mc
    else
        throw(ArgumentError("unknown p-value type."))
    end

    return p_val
end

"""
    fit(model::Function, xdata::AbstractArray, ydata::AbstractArray{AD.uwreal}, p0::AbstractArray; fit_type=:correlated_posdef, gaussian_priors=nothing, gof=false, kargs...) -> fit_result::FitResult

Perform a fit of the model function `model` to data `(xdata, ydata)`. The error is
automatically propagated to the fit parameters.

### Arguments
- `model(xdata::Vector, p::Vector)`: The model function to fit the data.
- `xdata::AbstractArray`: The independent variable data.
- `ydata::AbstractArray{AD.uwreal}`: The dependent variable data with uncertainties.
- `p0::AbstractArray`: The initial guess for the fit parameters.

### Keyword Arguments
- `fit_type`: The type of fit that is performed. Options are `:uncorrelated`, `:correlated`, or
  `:correlated_posdef`. With `fit_type=:correlated_posdef` (default) a positive definite
  approximation of the covariance matrix is used as the weight matrix. See the doc of the
  function `posdef_cov` for more informaiton.
- `gaussian_priors`: A `Dict{<:Integer, AD.uwreal}` or `Dict{<:Integer, <:AbstractVector}`
  which contains gaussian priors for the fit. The key must be the index of the corresponding
  parameter. In the first case the mean `μ` and error `σ` of the `AD.uwreal` object is used, 
  in the second case the `AbstractVector` is assumed to be of the form `[μ, σ]`.
  For each key `i` in the `Dict` the term `(p[i] - μ)²/σ²` is added to the chi-square
  function.
- `gof`: Whether to compute the goodness of fit (gof) metrics χ²_red and the p-value.
  Default is false.

### Returns
`fit_result`: A `FitResult` object containing the fit parameters and goodness of fit metrics.
It has the following fields:
- `param`::Vector{AD.uwreal}: The optimized parameters.
- `χ²_obs`: The observed `χ²` for the optimized parameters.
- `χ²_red`: The reduces `χ²` defined as `χ²_obs/χ²_exp` where `χ²_exp` is the expected
  chi-square value computed using `AD.chiexp`.
- `dof`: The numbers of degrees of freedom of the fit (`N_y - N_priors - N_params`)
- `p_value`: P-value of the fit computed as described in arXiv:2209.14188. See the doc of
  the function `p_value` for more information.
"""
function fit(model::Function, xdata::AbstractArray, ydata::AbstractArray{AD.uwreal},
             p0::AbstractArray; fit_type=:correlated_posdef,
             gaussian_priors=nothing, gof=false, kargs...)
    # Compute error
    err!.(ydata)

    # Ignore values that are not finite
    valid_indices = isfinite.(AD.value.(ydata))
    xdata_ = xdata[valid_indices]
    ydata_ = ydata[valid_indices]
    N_y = length(ydata_)

    # Get values and error of ydata
    ydata_value = AD.value.(ydata_)
    ydata_err = AD.err.(ydata_)

    # Define prior function
    if !isnothing(gaussian_priors)
        # Bring priors in correct format
        if typeof(gaussian_priors) <: Dict{<:Integer, AD.uwreal}
            gaussian_priors = Dict(k => [AD.value(a), AD.err(err!(a))]
                                   for (k, a) in gaussian_priors)
        elseif !(typeof(gaussian_priors) <: Dict{<:Integer, <:AbstractVector})
            throw(ArgumentError("gaussian_priors must be a Dict{<:Integer, AD.uwreal}, "*
                                "Dict{<:Integer, <:AbstractVector} or nothing."))
        end
        
        prior = (p) -> [(p[i] - μ)/σ for (i, (μ, σ)) in gaussian_priors]
        N_priors = length(gaussian_priors)
    else
        gaussian_priors = Dict{Int, Vector{Float64}}()
        prior = (p) -> []
        N_priors = 0
    end

    # Set weights and cost function
    if fit_type == :uncorrelated
        # Weights: Diag matrix of inverse of variances
        W = LA.diagm([(@. 1/ydata_err^2)...,
                      (1/σ^2 for (_, σ) in values(gaussian_priors))...])

        # Cost function
        u = @. 1/ydata_err
        cost = (p) -> [(u .* (model(xdata_, p) - ydata_value))..., prior(p)...]
    elseif fit_type in [:correlated, :correlated_posdef]
        # Weights
        if fit_type == :correlated
            # Covariance matrix
            W_ = inv(cov(ydata_))
        elseif fit_type == :correlated_posdef
            # "Positive definite" covariance matrix
            W_ = inv(posdef_cov(ydata_))
        end
        W_ = 0.5*(W_ + W_')

        # Extend weight matrix with variance of priors
        if N_priors != 0
            W_prior = LA.diagm([1/σ^2 for (_, σ) in values(gaussian_priors)])
            W = zeros(N_y + N_priors, N_y + N_priors)
            W[1:N_y, 1:N_y] = W_
            W[N_y+1:end, N_y+1:end] = W_prior
        else
            W = W_
        end

        # Cost function
        u = LA.cholesky(W_).U
        cost = (p) -> [(u * (model(xdata_, p) - ydata_value))..., prior(p)...]
    else
        throw(ArgumentError("unknown fit type. Use :uncorrelated, :correlated or "*
                            ":correlated_posdef."))
    end

    # Perform fit
    fit_result = LsqFit.lmfit(cost, p0, W)

    # Χ² function for propagating error to fit parameters
    function χ²(p, d)
        model_ydata_extended = [model(xdata_, p)...,
                                [p[i] for i in keys(gaussian_priors)]...]
        return (model_ydata_extended - d)' * W * (model_ydata_extended - d)
    end

    # Extend ydata with prior expectation values
    ydata_extended = [ydata_...,[AD.uwreal([μ, σ], "_prior$i")
                                 for (i, (μ, σ)) in gaussian_priors]...]

    # Observed χ² and number of degrees of freedom
    χ²_obs = χ²(fit_result.param, AD.value.(ydata_extended))
    dof = N_y + N_priors - length(p0)

    # Compute fit error, and reduced χ² and p-value if gof (goodness of fit) is true
    if gof
        fitp, χ²_exp = AD.fit_error(χ², fit_result.param, ydata_extended, parms.wpm,
                                    chi_exp=true)

        p_val = p_value(χ², ydata_extended, AD.value.(fitp), W, χ²_obs, dof,
                        fit_type=fit_type; kargs...)
        χ²_red = χ²_obs/χ²_exp
    else
        fitp = AD.fit_error(χ², fit_result.param, ydata_extended, parms.wpm, chi_exp=false)
        χ²_red = nothing
        p_val = nothing
    end

    fit_result = FitResult(fitp, χ²_obs, χ²_red, dof, p_val)
    return fit_result
end

"""
    fit_plateau(E_eff::AbstractVector{AD.uwreal}, plateau_range; guess=1.0, kargs...) -> fit_result::FitResult

Fit a constant function to the effective energy `E_eff` within the specified `plateau_range`.

### Arguments
- `E_eff::AbstractVector{AD.uwreal}`: The effective energy data.
- `plateau_range`: The range of indices within `E_eff` to fit the constant function to.
  It must be a `Vector` of length two of the form `[i_first, i_last]`.
- `guess=1.0`: Initial guess for the constant value.
- `kargs...`: Additional keyword arguments to be passed to the `fit` function.

### Returns
`fit_result::FitResult`: The result of the constant fit. See the doc of the function `fit`
for more information on its fields. 
"""
function fit_plateau(E_eff::AbstractVector{AD.uwreal}, plateau_range; guess=1.0, kargs...)
    # Fit to effective energy
    model(x, p) = @. p[1] + 0*x
    
    xdata = range(plateau_range...)
    ydata = E_eff[xdata.+1]
    p0 = [guess]

    return fit(model, xdata, ydata, p0; kargs...)
end