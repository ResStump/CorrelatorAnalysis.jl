"""
    read_uwreal(filename; return_info=false) -> a::Union{AD.uwreal, AbstractArray{AD.uwreal} [, info::String]

Read `AD.uwreal` objects from a BDIO file and return them as an array with the shape
specified in the file. If `length(a)==1` return it as a single object. The window
parameters are set in `parms.wpm` accordingly.
If `return_info=true`, additionally return the `info` string containing the creation
date and program versions.

## Note
The file is assumed to have the shape specified in the doc of the `write_uwreal`
function.
"""
function read_uwreal(filename; return_info=false)
    # Open BDIO file
    fb = BDIO.BDIO_open(string(filename), "r")

    # Read info, size of stored array, ensemble ids and window parameters
    info = nothing
    dims = nothing
    ids = nothing
    wpm = nothing
    while BDIO.BDIO_seek!(fb) && (info === nothing || dims === nothing)
        # Read info
        if BDIO.BDIO_get_uinfo(fb) == 0
            info = BDIO.BDIO_read_str(fb)
        end

        # Read size
        if BDIO.BDIO_get_uinfo(fb) == 1
            # Read number of dimensions
            ndims = Vector{Int64}(undef, 1)
            BDIO.BDIO_read(fb, ndims)

            # Read size of array
            dims = Vector{Int64}(undef, ndims[1])
            BDIO.BDIO_read(fb, dims)

            # Read ensemble ids and window parameters
            nids = Vector{Int64}(undef, 1)
            BDIO.BDIO_read(fb, nids)
            if nids[1] > 0
                ids = Vector{Int64}(undef, nids[1])
                BDIO.BDIO_read(fb, ids)
                wpm = Vector{Float64}(undef, 4*nids[1])
                BDIO.BDIO_read(fb, wpm)
            else
                ids = Int[]
                wpm = Float64[]
            end
        end
    end
    if info === nothing || dims === nothing
        throw(ArgumentError("the dimension are not stored in the right format."))
    end

    # Read uwreal objects
    len = prod(dims)
    idx = 0
    a = Vector{AD.uwreal}(undef, len)
    while BDIO.BDIO_seek!(fb) && idx < len
        if BDIO.BDIO_get_uinfo(fb) == 2
            idx += 1
            a[idx] = AD.read_uwreal(fb)
        end
    end
    if idx != len
        throw(ArgumentError("file does not contain the full AD.uwreal array."))
    end
    
    BDIO.BDIO_close!(fb)

    # Reshape array
    a = reshape(a, dims...)
    if length(a) == 1
        a = a[1]
    end

    # Set window parameters
    wpm = reshape(wpm, 4, :)
    for (i, id) in enumerate(ids)
        ensemble = AD.get_name_from_id(id, AD.wsg)
        parms.wpm[ensemble] = wpm[:, i]
    end

    if return_info
        return a, info
    else
        return a
    end
end

"""
    write_uwreal(filename, a::AD.uwreal, mode="d")
    write_uwreal(filename, a::AbstractArray{AD.uwreal}, mode="d")

Write the `AD.uwreal` or array of `AD.uwreal` `a` to a BDIO file. Specify how the file is
opend with `mode` (default is "d").

Additionally, write an info string containing the creation date and program version together
with the shape of the array and the relevant window parameters in `parms.wpm` to the file.

### Possible modes
- Write mode ("w"): The file is created. If the file exists an error is
  printed.
- Delete mode ("d"): The file is created. If the file exists it is
  overwritten.

### Format of file
The info string, the shape, the window parameters and the `AD.uwreal` are stored to
different records in the following way:
- Record 0: Info string as `BDIO_ASC_GENERIC`.
- Record 1: ndims(a), size(a)..., length(ids), ids..., wps... as `BDIO_BIN_GENERIC` where
    `ids::Vector{Int}` are the ensemble ids and `wps::Vector{Float64}` are the (flattened)
    window parameters relevent for `a`.
    (or set ndims(a) = 1 size(a) = 1 if it's just a single
  `AD.uwreal`).
- Record 2: The `AD.uwreal` are stored in consecutive records 2 using the function
  `AD.write_uwreal`.
"""
function write_uwreal(filename, a::AbstractArray{AD.uwreal}, mode="d")
    # Check that the mode is valid
    if mode ∉ ["w", "d"]
        throw(ArgumentError("mode must be 'w' or 'd'."))
    end

    # Create a BDIO file
    fb = BDIO.BDIO_open(string(filename), mode, "uwreal")

    # Write info string to record 0
    info = 
        "Date: $(Dates.now())\n"*
        "Julia version: $VERSION\n"*
        "$(@__MODULE__) version: $(pkgversion(@__MODULE__))\n"

    BDIO.BDIO_start_record!(fb, BDIO.BDIO_ASC_GENERIC, 0)
    BDIO.BDIO_write!(fb, info)

    BDIO.BDIO_write_hash!(fb)

    # Get all ensemble ids for `a` that are in parms.wpm and their window parameters 
    AD.unique_ids!.(a, [AD.wsg])
    ids = Int[]
    wpm = Float64[]
    for (i, v) in enumerate(a)
        for id in v.ids
            ensemble = AD.get_name_from_id(id, AD.wsg)
            if id ∉ ids && haskey(parms.wpm, ensemble)
                push!(ids, id)
                append!(wpm, parms.wpm[ensemble])
            end
        end

        # If v doesn't have an error, multiply by uwone
        if length(v.ids) == 0
            a[i] = err!(v*uwone)
        end
    end

    # Start writing to record 1
    BDIO.BDIO_start_record!(fb, BDIO.BDIO_BIN_GENERIC, 1)
    
    # Write size of array and number of ensembles
    BDIO.BDIO_write!(fb, Int[ndims(a), size(a)..., length(ids)])

    # Write the ensemble ids and then their window parameters
    if length(ids) > 0
        BDIO.BDIO_write!(fb, ids)
        BDIO.BDIO_write!(fb, wpm)
    end

    BDIO.BDIO_write_hash!(fb)

    # Write the uwreal objects in `a` to record 2
    for v in a
        AD.write_uwreal(v, fb, 2)
    end

    BDIO.BDIO_close!(fb)
end
write_uwreal(filename, a::AD.uwreal, mode="d") = write_uwreal(filename, [a], mode)

"""
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