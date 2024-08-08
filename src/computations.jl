function err!(a::AD.uwreal)
    try
        AD.uwerr(a, parms.wpm)
    catch e
        println("Failed to calculate error.")
        println(e)
    end

    return a
end

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

function cov(a::AbstractVector{AD.uwreal})
    # Compute covariance matrix using the set window parameters wpm
    # This changes the error in `a` therfore compute it again using err!
    cov = AD.cov(a, parms.wpm)

    err!(a)

    return cov
end

function posdef_cov(a::AbstractVector{AD.uwreal})
    # Generate window parameters with window 1 (i.e. no autocorrelation)
    wpm = deepcopy(parms.wpm)
    for k in keys(wpm)
        wpm[k] = [1.0, -1.0, -1.0, -1.0]
    end

    # Get correlation matrix without autocorrelation (changes error in `a`!)
    cov0 = AD.cov(a, wpm)

    # Compute error of `a` including autocorrelation
    err!.(a)

    # Compute offdiagonal correlator matrix entries by assuming the integrated
    # autocorrelation time τ_ij of entrie [i, j] is the geometric mean of the
    # autocorrelation time of the observables `a[i]` and `a[j]`, i.e. τ_ij = √(τ_i*τ_j).
    # This results in a positive definite correlation matrix.
    autocorr_correction = AD.err.(a) ./ .√LA.diag(cov0)
    
    return LA.diagm(autocorr_correction) * cov0 * LA.diagm(autocorr_correction)
end

function effective_mass(Cₜ::AbstractVector{AD.uwreal}, variant=:log; guess=1.0,
                        folded=false)
    # Compute error
    err!.(Cₜ)

    if variant in [:exp, :log]
        m_eff = log.(Cₜ./circshift(Cₜ, -1))
    elseif variant in [:cosh, :sinh]
        N_elements = length(Cₜ)
        if folded
            Nₜ = (N_elements - 1) * 2
        else
            Nₜ = N_elements
        end
        
        m_eff = Vector{AD.uwreal}(undef, N_elements)

        # Define fit models
        f_cosh(m, p) = cosh(m*(p[1] - Nₜ/2)) / cosh(m*(p[1]+1 - Nₜ/2)) - p[2]
        f_sinh(m, p) = sinh(m*(p[1] - Nₜ/2)) / sinh(m*(p[1]+1 - Nₜ/2)) - p[2]

        if variant == :cosh
            f = f_cosh
        else
            f = sinh
        end

        # Loop over all times `t` and compute effective mass
        for (iₜ, t) in enumerate(0:N_elements-1)
            p = [AD.uwreal(Float64(t)), Cₜ[iₜ]/Cₜ[mod1(iₜ+1, N_elements)]]
            # Try to find root, otherwise set value to NaN
            try
                m_eff[iₜ] = abs(AD.root_error(f, Float64(guess), p))
            catch e
                m_eff[iₜ] = AD.uwreal(NaN)
                println("For t = $t (iₜ=$iₜ): ", e)
            end
        end
    else
        throw(ArgumentError("unknown variant to compute effective mass."))
    end

    return m_eff
end

function fit(model, xdata::AbstractArray, ydata::AbstractArray{AD.uwreal},
             p0::AbstractArray; chi_exp=true, correlated_fit=true)
    # Compute error
    err!.(ydata)

    # Ignore values that are not finite
    valid_indices = isfinite.(AD.value.(ydata))
    xdata_ = xdata[valid_indices]
    ydata_ = ydata[valid_indices]

    # Set weights
    if correlated_fit
        W = inv(posdef_cov(ydata_))
        W = 0.5*(W + W')
    else
        # Diag matrix of inverse of variance
        W = LA.diagm(@. 1/AD.err(ydata_)^2)
    end

    fit_result = LsqFit.curve_fit(model, xdata_, AD.value.(ydata_), W, p0)

    # Χ² function for propagating error to fit parameters
    function Χ²(p, d)
        model_ydata = model(xdata_, p)
        return (d - model_ydata)' * W * (d - model_ydata)
    end

    if chi_exp
        fitp, cexp = AD.fit_error(Χ², fit_result.param, ydata_, parms.wpm, chi_exp=true)
        return fitp, cexp
    else
        fitp = AD.fit_error(Χ², fit_result.param, ydata_, parms.wpm, chi_exp=false)
        return fitp
    end
end

function fit_plateau(m_eff::AbstractVector{AD.uwreal}, plateau_range; guess=1.0,
                     chi_exp=false, kargs...)
    # Fit to effective mass
    model(x, p) = @. p[1] + 0*x
    
    xdata = range(plateau_range...)
    ydata = m_eff[xdata.+1]
    p0 = [guess]

    if chi_exp
        fitp, cexp = fit(model, xdata, ydata, p0; chi_exp=true, kargs...)
        m = fitp[1]
        return m, cexp
    else
        fitp = fit(model, xdata, ydata, p0; chi_exp=false, kargs...)
        m = fitp[1]
        return m
    end
end