@testset "IO test" begin
    # Generate MC ensemble
    rng = Random.MersenneTwister(15)
    N, N_mc = 10, 1000

    τ = 4
    μ = 1.1
    σ_rel = 0.01
    samples = Array{Float64}(undef, N, N_mc)
    for i in 1:N
        samples[i, :] = CA.markov_chain(rng, N_mc, μ, μ*σ_rel, τ)
    end

    # Create uwreal array
    uwarr = CA.uwreal_array(samples, "MC ensemble", :auto)
    CA.err!.(uwarr)

    # Write and read single AD.uwreal object
    CA.write_uwreal("test.bdio", uwarr[1])
    a_read = CA.err!(CA.read_uwreal("test.bdio"))
    @test AD.value(a_read) ≈ AD.value(uwarr[1])
    @test AD.err(a_read) ≈ AD.err(uwarr[1])

    # Write AD.uwreal array
    CA.write_uwreal("test.bdio", uwarr)

    # Change wpm entry to see if it is set correctly
    wpm_old = CA.parms.wpm
    CA.parms.wpm["MC ensemble"] = Float64[0, 0, 0, 0]

    # Read AD.uwreal array
    uwarr_read = CA.err!(CA.read_uwreal("test.bdio"))
    @test AD.value.(uwarr_read) ≈ AD.value.(uwarr)
    @test AD.err.(uwarr_read) ≈ AD.err.(uwarr)

    # Check if window parameters are set correctly
    @test all(AD.ensembles.(uwarr) .== AD.ensembles.(uwarr_read))
    @test wpm_old == CA.parms.wpm

    # Write and read AD.uwreal array with reshaped data
    shape = (2, 5)
    CA.write_uwreal("test.bdio", reshape(uwarr, shape))
    uwarr_read = CA.err!(CA.read_uwreal("test.bdio"))
    @test AD.value.(uwarr_read) ≈ AD.value.(reshape(uwarr, shape))
    @test AD.err.(uwarr_read) ≈ AD.err.(reshape(uwarr, shape))

    # Write and read NaN
    CA.write_uwreal("test.bdio", AD.uwreal(NaN))
    CA.err!(CA.read_uwreal("test.bdio"))

    rm("test.bdio")

    
    # Test if CA.export_samples propagats delta exactly under an affine transformation
    A = rand(N, N)
    b = rand(N)
    obs1 = A*uwarr .+ b

    exported_samples = stack(CA.export_samples.(obs1), dims=1)
    @test exported_samples ≈ A*samples .+ b

    # Test that there is no bias from propagating the delta
    f(x) = exp(x)
    obs2 = f(uwarr[1])
    @test sum(CA.export_samples(obs2))/N_mc ≈ f(uwarr[1].mean)
end