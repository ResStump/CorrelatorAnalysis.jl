PrecompileTools.@setup_workload begin

# Input of uwreal's
a = AD.uwreal(rand(1000), "White noise")
b = AD.uwreal([1.0, 0.1], "Var with error")
p = AD.cobs([1.0, 2.0], [1.0 0.1; 0.1 2.0], "Parameters")

PrecompileTools.@compile_workload begin
    # ADerrors function calls
    #########################

    for op in (:-, :sin, :cos, :log, :log10, :log2, :sqrt, :exp, :exp2, :exp10, :sinh, :cosh, :tanh)
        @eval c = $op(a)
    end
    c = 1.0 + b
    c = 1.0 - b
    c = 1.0 * b
    c = 1.0 / b
    c = a + 2.0
    c = a - 2.0
    c = a * 2.0
    c = a / 2.0
    c = a ^ 3
    c = a + b
    c = a - b
    c = a * b
    c = a / b

    # Error analysis
    AD.uwerr(c)
    AD.cov([c, a, b])
    AD.trcov([1.0 2.0 3.0;
              2.0 1.0 0.4;
              3.0 0.4 1.0], [c, a, b])
    redirect_stdout(devnull) do 
        AD.details(c)
    end

    # Error analysis of fit parameters
    npt = 12
    sig = zeros(npt, npt)
    dx  = zeros(npt)
    for i in 1:npt
        dx[i]    = 0.01*i
        sig[i,i] = dx[i]^2
        for j in i+1:npt
            sig[i,j] = 0.0001 - 0.000005*abs(i-j)
            sig[j,i] = 0.0001 - 0.000005*abs(i-j)
        end
    end
    y = [0.0802273592699947
         0.09150837606934502
         0.047923648239388834
         0.024851583326401416
         0.01635482325054799
         0.13115281737744588
         0.21013177679604178
         0.002143355151617357
         0.2292183950698425
         -0.05174734852593241
         0.1384913891139784
         -0.05211234898997283]
    dt = AD.cobs(y, sig, "Fit data")
    chisq(p, d) = sum( (d .- p[1]) .^ 2 ./ dx .^2 )
    xp = [sum(AD.value.(dt) ./ dx)/sum(1.0 ./ dx)]
    fitp, csqexp = AD.fit_error(chisq, xp, dt)
    AD.chiexp(chisq, xp, dt)

    
    # CorrelatorAnalysis function calls
    ###################################

    # Create random correlator
    rng = Random.MersenneTwister(154)
    Nₜ, N_cnfg = 16, 100

    τ = 4
    A_exact, m_exact = 1e-2, 0.1
    σ_rel = 0.01

    corr = Array{Float64}(undef, Nₜ, N_cnfg)
    for (iₜ, t) in enumerate(0:Nₜ-1)
        μ = A_exact*cosh(m_exact*(t - Nₜ/2))
        corr[iₜ, :] = markov_chain(rng, N_cnfg, μ, σ_rel*μ, τ)
    end

    # Initialize correlator and fold it
    Cₜ = uwreal_array(corr, "random", 1)
    Cₜ = fold_correlator(Cₜ)

    # Compute effective mass
    am_eff = effective_mass(Cₜ, :cosh, folded=true)

    # Constant fit to plateau
    plateau_range = [2, 6]
    fit_result = fit_plateau(am_eff, plateau_range)
    am = fit_result.param[1]
    err!(am)

    # Fit to correlator
    corr_model(t, p) = @. p[1]*cosh(p[2]*(t - Nₜ/2))
    fit_range = [plateau_range[1], plateau_range[2]+1]
    xdata = range(fit_range...)
    ydata = Cₜ[xdata.+1]

    p0 = [5e-3, 0.11]
    fit_result = fit(corr_model, xdata, ydata, p0)

    # Plot function calls
    plot_correlator(Cₜ)
    plot_model!(corr_model, [0, Nₜ÷2], AD.value.([A, am_fit]))
    plot_effective_mass(am_eff)
    plot_error_rectangle!(am, plateau_range)
end
end