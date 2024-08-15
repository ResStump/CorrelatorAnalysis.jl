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