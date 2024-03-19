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

function effective_mass(Cₜ::Vector{AD.uwreal}, variant=:log; guess=1.0, folded=false)
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

function fit(model, xdata::AbstractArray, ydata::Array{AD.uwreal}, p0::AbstractArray;
             chi_exp=true)
    # Compute error
    err!.(ydata)

    # Ignore values that are not finite
    valid_indices = isfinite.(AD.value.(ydata))
    xdata_ = xdata[valid_indices]
    ydata_ = ydata[valid_indices]
    W = @. 1/AD.err(ydata_)^2

    fit_result = LsqFit.curve_fit(model, xdata_, AD.value.(ydata_), W, p0)

    # Χ² function
    function Χ²(p, d)
        model_ydata = model(xdata_, p)
        return sum(@. (d - model_ydata)^2 * W)
    end

    if chi_exp
        fitp, cexp = AD.fit_error(Χ², fit_result.param, ydata_, parms.wpm, chi_exp=true)
        return fitp, cexp
    else
        fitp = AD.fit_error(Χ², fit_result.param, ydata_, parms.wpm, chi_exp=false)
        return fitp
    end
end

function fit_plateau(m_eff::Vector{AD.uwreal}, plateau_range; guess=1.0, chi_exp=false)
    # Fit to effective mass
    model(x, p) = @. p[1] + 0*x
    
    xdata = range(plateau_range...)
    ydata = m_eff[xdata.+1]
    p0 = [guess]

    if chi_exp
        fitp, cexp = fit(model, xdata, ydata, p0, chi_exp=true)
        m = fitp[1]
        return m, cexp
    else
        fitp = fit(model, xdata, ydata, p0, chi_exp=false)
        m = fitp[1]
        return m
    end
end