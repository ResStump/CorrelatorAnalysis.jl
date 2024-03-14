@testset "utils test" begin
    # Generate MC ensemble
    N, N_mc = 64, 1000
    eta  = randn(N_mc);
    x = Array{Float64}(undef, N, N_mc)

    # Random walk in [-1, 1]
    x[:, 1] .= 0.0
    for i in 2:1000
        accept = abs.(x[:, i]) .<= 1.0
        x[:, i] = @. x[:, i-1] + accept*eta[i]
    end

    # Initislize uwreal arrays
    mcid = "Random walk ensemble in [-1,1]"
    uwarr1 = [AD.uwreal(x[i, :], mcid) for i in 1:N]
    uwarr2 = CA.uwreal_array(x, mcid, :auto, :last)

    # Compute error
    CA.err!.(uwarr1)

    @test AD.value.(uwarr1) == AD.value.(uwarr2)
    @test AD.err.(uwarr1) == AD.err.(uwarr2)
end