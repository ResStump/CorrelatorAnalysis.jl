"""
    read_uwreal(file; return_info=false) -> A::Union{AD.uwreal, AbstractArray{AD.uwreal} [, info::String]

Read the `AD.uwreal` objects stored in `file` and return them as an array with the
specified shape. If `length(A)==1` return it as a single object.
If `return_info=true`, additionally return the `info` string containing the creaton
data and program versions.

## Note
The file is assumed to have the shape specified in the doc of the `write_uwreal`
function.
"""
function read_uwreal(file; return_info=false)
    # Open BDIO file
    fb = BDIO.BDIO_open(string(file), "r")

    # Read info and size of stored array
    info = nothing
    dims = nothing
    while BDIO.BDIO_seek!(fb) && (info === nothing || dims === nothing)
        # Read info
        if BDIO.BDIO_get_uinfo(fb) == 0
            info = BDIO.BDIO_read_str(fb)
        end

        # Read size
        if BDIO.BDIO_get_uinfo(fb) == 1
            # Read number of dimensions
            ndims = Vector{Int64}(undef, 1)
            BDIO.BDIO_read(fb, ndims, 1)

            # Read size of array
            dims = Vector{Int64}(undef, ndims[1])
            BDIO.BDIO_read(fb, dims)
        end
    end
    if info === nothing || dims === nothing
        throw(ArgumentError("the dimension are not stored in the right format."))
    end

    # Read uwreal objects
    len = prod(dims)
    idx = 0
    A = Vector{AD.uwreal}(undef, len)
    while BDIO.BDIO_seek!(fb) && idx < len
        if BDIO.BDIO_get_uinfo(fb) == 2
            idx += 1
            A[idx] = AD.read_uwreal(fb)
        end
    end
    if idx != len
        throw(ArgumentError("file does not contain the full AD.uwreal array."))
    end
    
    BDIO.BDIO_close!(fb)

    # Reshape array
    A = reshape(A, dims...)
    if length(A) == 1
        A = A[1]
    end

    if return_info
        return A, info
    else
        return A
    end
end

"""
    write_uwreal(file, a::AD.uwreal, mode="d")
    write_uwreal(file, a::AbstractArray{AD.uwreal}, mode="d")

Write the `AD.uwreal` or array of `AD.uwreal` `a` to the BDIO `file`. Specify how the file
is opend with `mode` (default is "d").

Additionally, write an info string containing the creation date and program version together
with the shape of the array to the file.

### Possible modes
- Write mode ("w"): The file is created. If the file exists an error is
  printed.
- Delete mode ("d"): The file is created. If the file exists it is
  overwritten.

### Format of file
The info string, the shape and the `AD.uwreal` are stored to different records in the
following way:
- Record 0: Info string as `BDIO_ASC_GENERIC`.
- Record 1: [ndims(a), size(a)...]  as `BDIO_BIN_INT64LE` (or [1, 1] if it's just a single
  `AD.uwreal`).
- Record 2: The `AD.uwreal` are stored in consecutive records 2 using the function
  `AD.write_uwreal`.
"""
function write_uwreal(file, a::AbstractArray{AD.uwreal}, mode="d")
    # Check that the mode is valid
    if mode ∉ ["w", "d"]
        throw(ArgumentError("mode must be 'w' or 'd'."))
    end

    # Create a BDIO file
    fb = BDIO.BDIO_open(string(file), mode, "uwreal")

    # Write info string to record 0
    info = 
        "Date: $(Dates.now())\n"*
        "Julia version: $VERSION\n"*
        "$(@__MODULE__) version: $(pkgversion(@__MODULE__))\n"

    BDIO.BDIO_start_record!(fb, BDIO.BDIO_ASC_GENERIC, 0)
    BDIO.BDIO_write!(fb, info)
    BDIO.BDIO_write_hash!(fb)

    # Write size of array to record 1
    BDIO.BDIO_start_record!(fb, BDIO.BDIO_BIN_INT64LE, 1)
    BDIO.BDIO_write!(fb, [ndims(a), size(a)...])
    BDIO.BDIO_write_hash!(fb)

    # Write the uwreal objects in `A` to record 2
    for v in a
        AD.write_uwreal(v, fb, 2)
    end

    BDIO.BDIO_close!(fb)
end
write_uwreal(file, a::AD.uwreal, mode="d") = write_uwreal(file, [a], mode)

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