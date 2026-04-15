@testset "utils test" begin
    # Generate MC ensemble
    rng = Random.MersenneTwister(12)
    N, N_mc = 64, 1000
    eta  = randn(rng, N_mc);
    x = Array{Float64}(undef, N, N_mc)

    # Random walk in [-1, 1]
    x[:, 1] .= 0.0
    for i in 2:1000
        accept = abs.(x[:, i]) .<= 1.0
        x[:, i] = @. x[:, i-1] + accept*eta[i]
    end

    # Initialize uwreal
    mcid = "Random walk ensemble in [-1,1]"
    a1 = AD.uwreal(x[1, :], mcid)
    a2 = CA.uwreal(x[1, :], mcid, :auto)

    # Compute error
    CA.err!(a1)

    @test a1.mean == a2.mean
    @test a1.err == a2.err
    @test a1 ≈ a2

    # Initialize uwreal arrays
    uwarr1 = [AD.uwreal(x[i, :], mcid) for i in 1:N]
    uwarr2 = CA.uwreal_array(x, mcid, :auto, :last)

    # Compute error
    CA.err!.(uwarr1)

    @test uwarr1 ≈ uwarr2

    # Propagate error using finite differences
    @test CA.derivedobs_fd.(exp, uwarr1) ≈ exp.(uwarr1)
    f = (a, x) -> a*exp(x)
    @test CA.derivedobs_fd.(f, 2, uwarr1) ≈ f.(2, uwarr1)

    # Test pencil_of_function for specific case
    C = [1 9 17; 2 10 18; 3 11 19; 4 12 20; 5 13 21; 6 14 22; 7 15 23; 8 16 24;;;
         25 33 41; 26 34 42; 27 35 43; 28 36 44; 29 37 45; 30 38 46; 31 39 47; 32 40 48;;;
         49 57 65; 50 58 66; 51 59 67; 52 60 68; 53 61 69; 54 62 70; 55 63 71; 56 64 72]
    C_pof = [1 2 9 11 17; 2 3 10 12 18; 3 4 11 13 19; 4 5 12 14 20;;;
             2 3 10 12 18; 3 4 11 13 19; 4 5 12 14 20; 5 6 13 15 21;;;
             25 26 33 35 41; 26 27 34 36 42; 27 28 35 37 43; 28 29 36 38 44;;;
             27 28 35 37 43; 28 29 36 38 44; 29 30 37 39 45; 30 31 38 40 46;;;
             49 50 57 59 65; 50 51 58 60 66; 51 52 59 61 67; 52 53 60 62 68]
    idx_pof = [1, 2]
    Δ = [1, 2]
    @test CA.pencil_of_function(C, idx_pof, Δ=Δ) == C_pof
end