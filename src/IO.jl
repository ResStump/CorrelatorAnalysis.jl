@doc raw"""
    read_uwreal(file) -> a

    Reads the `AD.uwreal` object from the `file`.
"""
function read_uwreal(file)
    # Open a BDIO file and move to first record
    fb = BDIO.BDIO_open(string(file), "r")
    BDIO.BDIO_seek!(fb)

    # Read the uwreal object from the file
    a = AD.read_uwreal(fb)

    BDIO.BDIO_close!(fb)

    return a
end

@doc raw"""
write_uwreal(a::AD.uwreal, file, mode="d")

Create a BDIO file and store the `AD.uwreal` object `a` in it.
"""
function write_uwreal(a::AD.uwreal, file, mode="d")
    # Create a BDIO file and store the uwreal object `a` in it
    fb = BDIO.BDIO_open(string(file), mode, "uwreal")
    AD.write_uwreal(a, fb, 8)
    BDIO.BDIO_close!(fb)
end

@doc raw"""
    export_samples(a::AD.uwreal) -> samples

Export the samples from the `AD.uwreal` object `a` which depends on exactly one ensemble.
The deviations from the mean of the primary observables are propagated linearly to `a` for
each MCMC sample and mean + deviation is returned.
"""
function export_samples(a::AD.uwreal)
    # Compute error
    err!(a)

    # Check that `a` depends on exactly one ensemble
    mcid = a.ids[1]
    idx = AD.wsg.map_ids[mcid]
    if length(a.ids) != 1 || AD.wsg.fluc[idx].nd == 1
        throw(DomainError("the uwreal object must depend on exactly one ensemble."))
    end

    # Propagate delta
    nd = AD.wsg.fluc[idx].nd
    new_delta  = zeros(nd)
    for i in eachindex(a.prop)
        if a.prop[i]
            new_delta .+= a.der[i]*AD.wsg.fluc[i].delta
        end
    end

    return a.mean .+ new_delta
end