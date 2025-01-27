@testset "Correlator 1" begin
    # Read openqxd correlator
    file_path = "data/correlator1.hdf5"
    corr = HDF5.h5read(file_path, "Correlator")
    Nₜ = size(corr)[1]

    # Lattice spacing
    a = AD.uwreal([0.0755, sqrt(0.0009^2 + 0.0007^2)], "a")/CA.ħc  # 1/MeV

    # Initialize correlator
    Cₜ = CA.uwreal_array(corr, "correlator1", 1)


    # Compute pion mass using non-folded correlator
    ###############################################

    # Compute effective mass
    am_eff = CA.effective_mass(Cₜ, :cosh)

    # Constant fit to plateau
    plateau_range = [13, 49]
    fit_result = CA.fit_plateau(am_eff, plateau_range, fit_type=:uncorrelated,
                                gof=false)
    am = fit_result.param[1]
    CA.err!(am)
    m = CA.err!(am/a)

    @test AD.value(am) ≈ 0.125145219671619
    @test AD.err(am) ≈ 0.0022876954542135463
    @test AD.value(m) ≈ 327.0798452545256
    @test AD.err(m) ≈ 7.755521136320138

    # Fit to correlator
    corr_model(t, p) = @. p[1]*cosh(p[2]*(t - Nₜ/2))
    fit_range = [plateau_range[1], plateau_range[2]+1]
    xdata = range(fit_range...)
    ydata = Cₜ[xdata.+1]

    p0 = [2e-3, 0.11]
    fit_result = CA.fit(corr_model, xdata, ydata, p0, fit_type=:uncorrelated, gof=false)
    (A, am_fit) = fit_result.param
    CA.err!(A)
    m_fit = CA.err!(am_fit/a)

    @test AD.value(m_fit) ≈ 325.299702527668
    @test AD.err(m_fit) ≈ 9.70021045013493
    @test AD.value(A) ≈ 0.002285894736669282
    @test AD.err(A) ≈ 0.00022145561505033227


    # Compute pion mass using folded correlator
    ###########################################

    Cₜ_folded = CA.fold_correlator(Cₜ)

    # Compute effective mass
    am_eff_folded = CA.effective_mass(Cₜ_folded, :cosh, folded=true)

    # Constant fit to plateau
    plateau_range_folded = [12, 31]
    fit_result = CA.fit_plateau(am_eff_folded, plateau_range_folded, fit_type=:uncorrelated,
                                gof=false)
    am_folded = fit_result.param[1]
    CA.err!(am_folded)
    m_folded = CA.err!(am_folded/a)

    @test AD.value(am_folded) ≈ 0.12507519414369922
    @test AD.err(am_folded) ≈ 0.002505118909142919
    @test AD.value(m_folded) ≈ 326.8968263673802
    @test AD.err(m_folded) ≈ 8.199945577057564

    fit_range_folded = [plateau_range_folded[1], plateau_range_folded[2]+1]
    xdata = range(fit_range_folded...)
    ydata = Cₜ_folded[xdata.+1]

    p0 = [5e-3, 0.11]
    fit_result = CA.fit(corr_model, xdata, ydata, p0, fit_type=:uncorrelated, gof=false)
    (A_folded, am_folded_fit) = fit_result.param
    CA.err!(A_folded)
    m_folded_fit = CA.err!(am_folded_fit/a)

    @test AD.value(m_folded_fit) ≈ 327.1734476838068
    @test AD.err(m_folded_fit) ≈ 9.119123879443574
    @test AD.value(A_folded) ≈ 0.002300828559365982
    @test AD.err(A_folded) ≈ 0.00023094851116168354
end

@testset "Correlator 2" begin
    # Read openqxd correlator
    file_path = "data/correlator2.hdf5"
    corr = HDF5.h5read(file_path, "Correlator")
    Nₜ, N_src, N_cnfg = size(corr)

    # Lattice spacing
    a = AD.uwreal([0.0762, 0.0], "a")/CA.ħc  # 1/MeV

    # Initialize correlator
    corr_src_mean = sum(real(corr), dims=2)[:, 1, :]/N_src
    Cₜ = CA.uwreal_array(corr_src_mean, "correlator2", :auto)


    # Compute pion mass using non-folded correlator
    ###############################################

    # Compute effective mass
    am_eff = CA.effective_mass(Cₜ, :cosh)

    # Constant fit to plateau
    plateau_range = [11, 51]
    fit_result = CA.fit_plateau(am_eff, plateau_range, fit_type=:uncorrelated,
                                gof=false)
    am = fit_result.param[1]
    m = CA.err!(am/a)

    @test AD.value(m) ≈ 657.3138400277402
    @test AD.err(m) ≈ 1.5780809793144068

    # Fit to correlator
    corr_model(t, p) = @. p[1]*cosh(p[2]*(t - Nₜ/2))
    fit_range = [plateau_range[1], plateau_range[2]+1]
    xdata = range(fit_range...)
    ydata = Cₜ[xdata.+1]

    p0 = [7e-3, 0.25]
    fit_result = CA.fit(corr_model, xdata, ydata, p0, fit_type=:uncorrelated, gof=false)
    (A, am_fit) = fit_result.param
    CA.err!(A)
    m_fit = CA.err!(am_fit/a)

    @test AD.value(m_fit) ≈ 656.5952337578769
    @test AD.err(m_fit) ≈ 1.6734785131400542
    @test AD.value(A) ≈ 0.006512432451529093
    @test AD.err(A) ≈ 0.00010871594531586655


    # Compute pion mass using folded correlator
    ###########################################

    Cₜ_folded = CA.fold_correlator(Cₜ)

    # Compute effective mass
    am_eff_folded = CA.effective_mass(Cₜ_folded, :cosh, folded=true)

    # Constant fit to plateau
    plateau_range_folded = [11, 28]
    fit_result = CA.fit_plateau(am_eff_folded, plateau_range_folded, fit_type=:uncorrelated,
                                gof=false)
    am_folded = fit_result.param[1]
    m_folded = CA.err!(am_folded/a)


    @test AD.value(m_folded) ≈ 657.6202561282212
    @test AD.err(m_folded) ≈ 1.5212855009840285

    fit_range_folded = [plateau_range_folded[1], plateau_range_folded[2]+1]
    xdata = range(fit_range_folded...)
    ydata = Cₜ_folded[xdata.+1]

    p0 = [7e-3, 0.25]
    fit_result = CA.fit(corr_model, xdata, ydata, p0, fit_type=:uncorrelated,
                                gof=false)
    (A_folded, am_folded_fit) = fit_result.param
    CA.err!(A_folded)
    m_folded_fit = CA.err!(am_folded_fit/a)

    @test AD.value(m_folded_fit) ≈ 656.8686619785077
    @test AD.err(m_folded_fit) ≈ 1.5926570085408367
    @test AD.value(A_folded) ≈ 0.006503071705487689
    @test AD.err(A_folded) ≈ 0.0001060639274430725
end

@testset "p-value" begin
    # Creat random correlator
    # (with large number of cnfgs to get positive definite covariance matrix)
    rng = Random.MersenneTwister(154)
    Nₜ, N_cnfg = 16, 10000

    τ = 4
    A_exact, m_exact = 1e-2, 0.1
    σ_rel = 0.005

    corr = Array{Float64}(undef, Nₜ, N_cnfg)
    for (iₜ, t) in enumerate(0:Nₜ-1)
        μ = A_exact*cosh(m_exact*(t - Nₜ/2))
        corr[iₜ, :] = CA.markov_chain(rng, N_cnfg, μ, σ_rel*μ, τ)
    end

    mcid = "random_10000"
    Cₜ = CA.uwreal_array(corr, mcid, :auto)

    # Fold correlator
    Cₜ_folded = CA.fold_correlator(Cₜ)

    # Effective mass of folded correlator
    am_eff_folded = CA.effective_mass(Cₜ_folded, :cosh, folded=true)
    CA.err!.(am_eff_folded)

    # Perform correlated fit and compute p-value once with exact formula and once with 
    # general formula using MC
    plateau_range_folded = [2, 6]
    p_value_exact = CA.fit_plateau(am_eff_folded, plateau_range_folded,
                                   fit_type=:correlated, gof=true,
                                   p_value_type=:correlated).p_value
    p_value_general = CA.fit_plateau(am_eff_folded, plateau_range_folded,
                                     fit_type=:correlated, gof=true,
                                     p_value_type=:general, N_mc=10^6, rng=rng).p_value

    # Compare p-values to 3 decimal palces
    @test abs(p_value_exact - p_value_general) < 1e-3
end

@testset "GEVP" begin
    # Read data and initialize correlator matrix
    file_path = "data/correlator_matrix.hdf5"
    file = HDF5.h5open(file_path)
    corr_matrix = read(file["Correlator Matrix"])
    close(file)
    Nₜ, N_op, _, N_cnfg = size(corr_matrix)
    Cₜ = CA.uwreal_array(corr_matrix, "B450r000", :auto)

    # Compute effective energy using t₀ = ceil(t/2) method
    plateau_range = [15, 24]
    aE_eff = CA.GEVP(Cₜ, :ceil_t_half)

    # Extract lowest 3 energy levels
    fit_result = CA.fit_plateau(aE_eff[1], plateau_range, fit_type=:correlated_posdef)
    aE0 = CA.err!(fit_result.param[1])
    fit_result = CA.fit_plateau(aE_eff[2], plateau_range, fit_type=:correlated_posdef)
    aE1 = CA.err!(fit_result.param[1])
    fit_result = CA.fit_plateau(aE_eff[3], plateau_range, fit_type=:correlated_posdef)
    aE2 = CA.err!(fit_result.param[1])

    @test AD.value(aE0) ≈ 2.1524622389214803
    @test AD.err(aE0) ≈ 0.0008896215919601967
    @test AD.value(aE1) ≈ 2.17874053787349
    @test AD.err(aE1) ≈ 0.0009186633907838645
    @test AD.value(aE2) ≈ 2.1883420746340887
    @test AD.err(aE2) ≈ 0.0008906508771061862

    # Compute effective energy using t₀ = const method
    t₀ = 12
    plateau_range = [15, 24]
    aE_eff = CA.GEVP(Cₜ, t₀)

    # Extract lowest 3 energy levels
    fit_result = CA.fit_plateau(aE_eff[1], plateau_range, fit_type=:correlated_posdef)
    aE0 = CA.err!(fit_result.param[1])
    fit_result = CA.fit_plateau(aE_eff[2], plateau_range, fit_type=:correlated_posdef)
    aE1 = CA.err!(fit_result.param[1])
    fit_result = CA.fit_plateau(aE_eff[3], plateau_range, fit_type=:correlated_posdef)
    aE2 = CA.err!(fit_result.param[1])

    @test AD.value(aE0) ≈ 2.152436415887633
    @test AD.err(aE0) ≈ 0.0008891187819048125
    @test AD.value(aE1) ≈ 2.1787670486036763
    @test AD.err(aE1) ≈ 0.0009180396787471865
    @test AD.value(aE2) ≈ 2.188342896037995
    @test AD.err(aE2) ≈ 0.0008907154702271366
end
