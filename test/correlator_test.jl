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

    # Compute effective energy
    am_eff = CA.effective_energy(Cₜ, :cosh)

    # Constant fit to plateau
    plateau_range = [13, 49]
    am = CA.fit_plateau(am_eff, plateau_range)
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
    (A, am_fit), cexp = CA.fit(corr_model, xdata, ydata, p0)
    CA.err!(A)
    m_fit = CA.err!(am_fit/a)

    @test AD.value(m_fit) ≈ 325.299702527668
    @test AD.err(m_fit) ≈ 9.70021045013493
    @test AD.value(A) ≈ 0.002285894736669282
    @test AD.err(A) ≈ 0.00022145561505033227


    # Compute pion mass using folded correlator
    ###########################################

    Cₜ_folded = CA.fold_correlator(Cₜ)

    # Compute effective energy
    am_eff_folded = CA.effective_energy(Cₜ_folded, :cosh, folded=true)

    # Constant fit to plateau
    plateau_range_folded = [12, 31]
    am_folded = CA.fit_plateau(am_eff_folded, plateau_range_folded)
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
    (A_folded, am_folded_fit), cexp = CA.fit(corr_model, xdata, ydata, p0)
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

    # Compute effective energy
    am_eff = CA.effective_energy(Cₜ, :cosh)

    # Constant fit to plateau
    plateau_range = [11, 51]
    am = CA.fit_plateau(am_eff, plateau_range)
    m = CA.err!(am/a)

    @test AD.value(m) ≈ 657.3138400277402
    @test AD.err(m) ≈ 1.5780809793144068

    # Fit to correlator
    corr_model(t, p) = @. p[1]*cosh(p[2]*(t - Nₜ/2))
    fit_range = [plateau_range[1], plateau_range[2]+1]
    xdata = range(fit_range...)
    ydata = Cₜ[xdata.+1]

    p0 = [7e-3, 0.25]
    (A, am_fit), cexp = CA.fit(corr_model, xdata, ydata, p0)
    CA.err!(A)
    m_fit = CA.err!(am_fit/a)

    @test AD.value(m_fit) ≈ 656.5952337578769
    @test AD.err(m_fit) ≈ 1.6734785131400542
    @test AD.value(A) ≈ 0.006512432451529093
    @test AD.err(A) ≈ 0.00010871594531586655


    # Compute pion mass using folded correlator
    ###########################################

    Cₜ_folded = CA.fold_correlator(Cₜ)

    # Compute effective energy
    am_eff_folded = CA.effective_energy(Cₜ_folded, :cosh, folded=true)

    # Constant fit to plateau
    plateau_range_folded = [11, 28]
    am_folded = CA.fit_plateau(am_eff_folded, plateau_range_folded)
    m_folded = CA.err!(am_folded/a)


    @test AD.value(m_folded) ≈ 657.6202561282212
    @test AD.err(m_folded) ≈ 1.5212855009840285

    fit_range_folded = [plateau_range_folded[1], plateau_range_folded[2]+1]
    xdata = range(fit_range_folded...)
    ydata = Cₜ_folded[xdata.+1]

    p0 = [7e-3, 0.25]
    (A_folded, am_folded_fit), cexp = CA.fit(corr_model, xdata, ydata, p0)
    CA.err!(A_folded)
    m_folded_fit = CA.err!(am_folded_fit/a)

    @test AD.value(m_folded_fit) ≈ 656.8686619785077
    @test AD.err(m_folded_fit) ≈ 1.5926570085408367
    @test AD.value(A_folded) ≈ 0.006503071705487689
    @test AD.err(A_folded) ≈ 0.0001060639274430725
end
