function plot_correlator!(p::Plots.Plot, Cₜ::Vector{AD.uwreal}; kargs...)
    # Compute error
    err!.(Cₜ)

    # Extract correlator values and set entries <= 0 to NaN
    Cₜ_values = AD.value.(Cₜ)
    Cₜ_values[Cₜ_values.<=0] .= NaN

    # Extract correlator errors
    Cₜ_err = AD.err.(Cₜ)

    # Get correlator label from ensemble names
    names = String[]
    for mcid in keys(parms.wpm)
        if AD.get_id_from_name(mcid) in Cₜ[1].ids
            push!(names, mcid)
        end
    end
    corr_label = join(names, ", ")

    Nₜ = length(Cₜ)
    Plots.plot!(p, yscale=:log, xlabel=L"t/a", ylabel=L"C(t)", legend=true,
                minorticks=true)
    Plots.scatter!(p, 0:Nₜ-1, Cₜ_values, yerror=Cₜ_err, label=corr_label,
                   markersize=2.5; kargs...)

    return p
end
plot_correlator(Cₜ::Vector{AD.uwreal}; kargs...) = 
    plot_correlator!(Plots.plot(), Cₜ; kargs...)
plot_correlator!(Cₜ::Vector{AD.uwreal}; kargs...) = 
    plot_correlator!(Plots.plot!(), Cₜ; kargs...)

function plot_autocorrelation(obs::AD.uwreal, mcid, N_cnfg_max=nothing; kargs...)
    # Compute error
    err!(obs)

    # Extract choosen window and normalized autocorrelation from `obs`
    iw = AD.window(obs, mcid)
    r = AD.rho(obs, mcid)
    dr = AD.drho(obs, mcid)

    # Specify how far the autocorrelation should be plotted
    if N_cnfg_max === nothing
        N_cnfg_max = 2iw
    end

    # Plot normalized autocorrelation
    p = Plots.plot(xlabel=L"i_\mathrm{cnfg}", ylabel=L"\rho", minorticks=true)
    Plots.scatter!(0:N_cnfg_max-1, r[1:N_cnfg_max], yerror=dr[1:N_cnfg_max],
                   label="Normalized autocorrelation"; kargs...)

    # Plot choosen window
    if parms.wpm[mcid][2] > 0
        # If S is specified, add it to the label
        S = parms.wpm[mcid][2]
        Plots.vline!([iw-1], label="Chosen window (S = $S)")
    else
        Plots.vline!([iw-1], label="Chosen window")
    end

    # Horizontal line at zero
    Plots.hline!([0], color=:black, label=nothing)

    return p
end

function plot_effective_mass!(p::Plots.Plot, m_eff::Vector{AD.uwreal}; unit="lattice",
                              kargs...)
    # Compute error
    err!.(m_eff)

    # String that contains unit (empty if lattice units are used)
    if unit == "lattice"
        unit_str = ""
    else
        unit_str = "[$unit]"
    end

    # Plot effective mass
    Nₜ = length(m_eff)
    Plots.plot!(p, xlabel=L"t/a", ylabel=L"$m_\mathrm{eff}(t)$ %$unit_str", minorticks=true)
    Plots.scatter!(p, 0:Nₜ-1, AD.value.(m_eff), yerror=AD.err.(m_eff),
                   label="Effective mass", markersize=2.5; kargs...)

    return p
end
plot_effective_mass(m_eff::Vector{AD.uwreal}; kargs...) = 
    plot_effective_mass!(Plots.plot(), m_eff; kargs...)
plot_effective_mass!(m_eff::Vector{AD.uwreal}; kargs...) = 
    plot_effective_mass!(Plots.plot!(), m_eff; kargs...)

# Helper function
rectangle(x1, x2, y, h) = Plots.Shape([x1, x2, x2, x1], y .+ [-h,-h,h,h])

function plot_error_rectangle!(p::Plots.Plot, m::AD.uwreal, plateau_range; kargs...)
    # Compute error
    err!(m)
    
    # Plot value of mass as horizontal line
    Plots.plot!(p, plateau_range, ones(2)*AD.value(m), color=:red, z_order=:front;
                label="Fit result", kargs...)
    # Plot red rectangular area
    Plots.plot!(p, rectangle(plateau_range[1], plateau_range[2], AD.value(m), AD.err(m)),
                color=:red, opacity=0.4, z_order=:back, label=nothing)

    return p
end
plot_error_rectangle(m::AD.uwreal, plateau_range; kargs...) = 
    plot_error_rectangle!(Plots.plot(), m, plateau_range; kargs...)
plot_error_rectangle!(m::AD.uwreal, plateau_range; kargs...) = 
    plot_error_rectangle!(Plots.plot!(), m, plateau_range; kargs...)

function plot_model!(p::Plots.Plot, model, xdata_range::AbstractVector,
                     parms::AbstractArray; kargs...)
    xdata = LinRange(xdata_range..., 100)
    
    Plots.plot!(p, xdata, model(xdata, parms), linewidth=2, label="Fit result"; kargs...)

    return p
end
plot_model(model, xdata::AbstractVector, parms::AbstractArray; kargs...) = 
    plot_model!(Plots.plot(), model, xdata, parms; kargs...)
plot_model!(model, xdata::AbstractVector, parms::AbstractArray; kargs...) = 
    plot_model!(Plots.plot!(), model, xdata, parms; kargs...)
