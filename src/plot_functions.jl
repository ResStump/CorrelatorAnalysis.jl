function plot_correlator!(p::Plots.Plot, Cₜ::AbstractVector{AD.uwreal}; t_shift=0, kargs...)
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
    t_arr = (0:Nₜ-1) .+ t_shift
    Plots.plot!(p, yscale=:log10, xlabel=L"t/a", ylabel=L"C(t)", legend=true,
                minorticks=true)
    Plots.scatter!(p, t_arr, Cₜ_values, yerror=Cₜ_err, label=corr_label,
                   markersize=2.5; kargs...)

    return p
end
plot_correlator(Cₜ::AbstractVector{AD.uwreal}; kargs...) = 
    plot_correlator!(Plots.plot(), Cₜ; kargs...)
plot_correlator!(Cₜ::AbstractVector{AD.uwreal}; kargs...) = 
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

function plot_effective_mass!(p::Plots.Plot, m_eff::AbstractVector{AD.uwreal};
                              unit="lattice", t_shift=0, kargs...)
    # Compute error
    err!.(m_eff)

    # String that contains unit (empty if lattice units are used)
    if unit == "lattice"
        ylabel = L"am_\mathrm{eff}(t)"
    else
        ylabel = L"$m_\mathrm{eff}(t)$ [%$unit]"
    end

    # Plot effective mass
    Nₜ = length(m_eff)
    t_arr = (0:Nₜ-1) .+ t_shift
    Plots.plot!(p, xlabel=L"t/a", ylabel=ylabel, minorticks=true)
    Plots.scatter!(p, t_arr, AD.value.(m_eff), yerror=AD.err.(m_eff),
                   label="Effective mass", markersize=2.5; kargs...)

    return p
end
plot_effective_mass(m_eff::AbstractVector{AD.uwreal}; kargs...) = 
    plot_effective_mass!(Plots.plot(), m_eff; kargs...)
plot_effective_mass!(m_eff::AbstractVector{AD.uwreal}; kargs...) = 
    plot_effective_mass!(Plots.plot!(), m_eff; kargs...)

# Helper function
rectangle(x1, x2, y, h) = Plots.Shape([x1, x2, x2, x1], y .+ [-h,-h,h,h])

function plot_error_rectangle!(p::Plots.Plot, m::AD.uwreal, plateau_range;
                               fill_kargs=Dict(), color=:red, kargs...)
    # Compute error
    err!(m)
    
    # Plot value of mass as horizontal line
    Plots.plot!(p, plateau_range, ones(2)*m.mean; color=color, z_order=:front,
                label="Fit result", kargs...)
    # Plot red rectangular area
    Plots.plot!(p, rectangle(plateau_range[1], plateau_range[2], m.mean, m.err);
                color=color, opacity=0.3, lineopacity=0, z_order=:back, label=nothing,
                fill_kargs...)

    return p
end
plot_error_rectangle(m::AD.uwreal, plateau_range; kargs...) = 
    plot_error_rectangle!(Plots.plot(), m, plateau_range; kargs...)
plot_error_rectangle!(m::AD.uwreal, plateau_range; kargs...) = 
    plot_error_rectangle!(Plots.plot!(), m, plateau_range; kargs...)

function plot_herrorline!(p::Plots.Plot, m::AD.uwreal; color=:red, fill_kargs=Dict(),
                         kargs...)
    # Compute error
    err!(m)

    Plots.hline!(p, [m.mean]; color=color, kargs...)
    Plots.hspan!(p, m.mean .+ [-m.err, m.err]; color=color, opacity=0.3, lineopacity=0,
                 z_order=:back, label=nothing, fill_kargs...)

    return p
end
plot_herrorline(m::AD.uwreal; kargs...) = 
    plot_herrorline!(Plots.plot(), m; kargs...)
plot_herrorline!(m::AD.uwreal; kargs...) =
    plot_herrorline!(Plots.plot!(), m; kargs...)

function plot_model!(p::Plots.Plot, model::Function, xdata_range::AbstractVector,
                     parms::AbstractArray; kargs...)    
    Plots.plot!(p, xdata -> model(xdata, parms), linewidth=2, label="Fit result"; kargs...)

    return p
end
plot_model(model, xdata::AbstractVector, parms::AbstractArray; kargs...) = 
    plot_model!(Plots.plot(), model, xdata, parms; kargs...)
plot_model!(model, xdata::AbstractVector, parms::AbstractArray; kargs...) = 
    plot_model!(Plots.plot!(), model, xdata, parms; kargs...)
