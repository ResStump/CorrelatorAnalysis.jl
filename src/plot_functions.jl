function plot_correlator!(ax::CM.Axis, Cₜ::AbstractVector{AD.uwreal}; t_shift=0,
                          marker=:circle, kargs...)
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

    # Swich to log scale
    limits = ax.limits[]
    CM.limits!(ax, nothing, nothing, eps(1.0), nothing) # Set y-axis limits to avoid log(0)
    ax.yscale = CM.log10
    if limits == (nothing, nothing)
        CM.limits!(ax, nothing, nothing, nothing, nothing)
    else
        CM.limits!(ax, limits...)
    end

    # Set axis labels
    ax.xlabel = L"t/a"
    ax.ylabel = L"C(t)"
    
    # Plot data
    scatter_plot = CM.scatter!(ax, t_arr, Cₜ_values; label=corr_label, marker=marker,
                               kargs...)
    errorbars_plot = CM.errorbars!(ax, t_arr, Cₜ_values, Cₜ_err; label=corr_label, kargs...)
    CM.axislegend(ax, merge=true)

    return ax, scatter_plot, errorbars_plot
end
plot_correlator!(Cₜ::AbstractVector{AD.uwreal}; kargs...) = 
    plot_correlator!(CM.current_axis(), Cₜ; kargs...)
plot_correlator(Cₜ::AbstractVector{AD.uwreal}; kargs...) = begin
    f, _, _ = CM.lines(Float64[])
    plot_correlator!(Cₜ; kargs...)
    return f
end

function plot_autocorrelation!(obs::AD.uwreal, mcid; N_cnfg_max=nothing, marker=:circle,
                               kargs...)
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

    # Set axis labels
    ax = CM.current_axis()
    ax.xlabel=L"i_\mathrm{cnfg}"
    ax.ylabel=L"\rho"

    # Plot choosen window
    if parms.wpm[mcid][2] > 0
        # If S is specified, add it to the label
        S = parms.wpm[mcid][2]
        CM.vlines!(ax, [iw-1], label="Chosen window (S = $S)", color=:red)
    else
        CM.vlines!(ax, [iw-1], label="Chosen window", color=:red)
    end

    # Plot horizontal line at zero
    CM.hlines!(ax, [0], color=:black)

    # Plot normalized autocorrelation
    scatter_plot = CM.scatter!(ax, 0:N_cnfg_max-1, r[1:N_cnfg_max],
                               label="Normalized autocorrelation"; marker=marker, kargs...)
    errorbars_plot = CM.errorbars!(ax, 0:N_cnfg_max-1, r[1:N_cnfg_max], dr[1:N_cnfg_max],
                                   label="Normalized autocorrelation"; kargs...)
                                   
    CM.axislegend(ax, merge=true)

    return ax, scatter_plot, errorbars_plot
end
plot_autocorrelation(obs::AD.uwreal, mcid; N_cnfg_max=nothing, kargs...) = begin
    f, _, _ = CM.lines(Float64[])
    plot_autocorrelation!(obs::AD.uwreal, mcid; N_cnfg_max=N_cnfg_max, kargs...)
    return f
end

function plot_effective_energy!(ax::CM.Axis, E_eff::AbstractVector{AD.uwreal};
                                unit="lattice", t_shift=0, marker=:circle, kargs...)
    # Compute error
    err!.(E_eff)

    # String that contains unit (empty if lattice units are used)
    if unit == "lattice"
    ylabel = L"aE_\mathrm{eff}(t)"
    else
    ylabel = L"$E_\mathrm{eff}(t)$ [%$unit]"
    end

    # Plot effective energy
    Nₜ = length(E_eff)
    t_arr = (0:Nₜ-1) .+ t_shift

    # Set axis attributes
    ax.xlabel=L"t/a"
    ax.ylabel=ylabel

    # Plot data
    scatter_plot = CM.scatter!(ax, t_arr, AD.value.(E_eff),
                               label="Effective energy"; marker=marker, kargs...)
    errorbars_plot = CM.errorbars!(ax, t_arr, AD.value.(E_eff), AD.err.(E_eff), 
                                   label="Effective energy"; kargs...)
    CM.axislegend(ax, merge=true)

    return ax, scatter_plot, errorbars_plot
end
plot_effective_energy!(E_eff::AbstractVector{AD.uwreal}; kargs...) = 
    plot_effective_energy!(CM.current_axis(), E_eff; kargs...)
plot_effective_energy(E_eff::AbstractVector{AD.uwreal}; kargs...) = begin
    f, _, _ = CM.lines(Float64[])
    plot_effective_energy!(E_eff; kargs...)
    return f
end

function plot_error_rectangle!(ax::CM.Axis, E::AD.uwreal, plateau_range; color=:red,
                               kargs...)
    # Compute error
    err!(E)

    lines_plot = CM.lines!(ax, plateau_range, [E.mean], color=color; kargs...)
    errorband_plot = CM.band!(ax, plateau_range, E.mean-E.err, E.mean+E.err,
                              color=(color, 0.3); kargs...)
    CM.axislegend(ax, merge=true)

    # Bring error band to back
    CM.translate!(errorband_plot, 0, 0, -10)

    return ax, lines_plot, errorband_plot
end
plot_error_rectangle!(E::AD.uwreal, plateau_range; kargs...) = 
    plot_error_rectangle!(CM.current_axis(), E, plateau_range; kargs...)

function plot_herrorline!(ax::CM.Axis, E::AD.uwreal; color=:red, kargs...)
    # Compute error
    err!(E)

    hlines_plot = CM.hlines!(ax, [E.mean]; color=color, kargs...)
    hspan_plot = CM.hspan!(ax, [E.mean-E.err], [E.mean+E.err]; color=(color, 0.3), kargs...)
    CM.axislegend(ax, merge=true)

    # Bring error band to back
    CM.translate!(hspan_plot, 0, 0, -10)

    return ax, hlines_plot, hspan_plot
end
plot_herrorline!(E::AD.uwreal; kargs...) =
    plot_herrorline!(CM.current_axis(), E; kargs...)

#= function plot_model!(p::Plots.Plot, model::Function, xdata_range::AbstractVector,
                     parms::AbstractArray; kargs...)    
    Plots.plot!(p, xdata -> model(xdata, parms), linewidth=2, label="Fit result"; kargs...)

    return p
end
plot_model(model, xdata::AbstractVector, parms::AbstractArray; kargs...) = 
    plot_model!(Plots.plot(), model, xdata, parms; kargs...)
plot_model!(model, xdata::AbstractVector, parms::AbstractArray; kargs...) = 
    plot_model!(Plots.plot!(), model, xdata, parms; kargs...) =#

function plot_model!(ax::CM.Axis, model::Function, xdata_range::AbstractVector,
                     parms::AbstractArray; kargs...)    
    lines_plot = CM.lines!(ax, xdata_range[1]..xdata_range[end], xdata -> model(xdata, parms),
                           label="Fit result"; kargs...)
    CM.axislegend(ax, merge=true)

    return ax, lines_plot
end
plot_model!(model, xdata::AbstractVector, parms::AbstractArray; kargs...) = 
    plot_model!(CM.current_axis(), model, xdata, parms; kargs...)
plot_model(model, xdata::AbstractVector, parms::AbstractArray; kargs...) = begin
    f, _, _ = CM.lines(Float64[])
    plot_model!(model, xdata, parms; kargs...)
    return f
end

function plot_overlaps!(f::CM.Figure, Z_in::AbstractMatrix; n_range=:all,
                        colors=Makie.wong_colors())
    # Number of operators
    N_op = size(Z_in, 1)

    # Energy level range for plot
    if n_range == :all
        n_range = 1:size(Z_in, 2)
    end
    n_levels = length(n_range)

    # Arrange data such that each column corresponds the levels of one operator
    # which are stacked on top of each other
    bar_heights = transpose(Z_in[:, n_range])
    bar_xpos = stack([1:N_op for n in n_range], dims=1)
    groups = stack([n_range for i in 1:N_op], dims=2) # groups with same color

    # Flatten data
    bar_heights, bar_xpos, groups = vec(bar_heights), vec(bar_xpos), vec(groups)

    # Create Axis
    ax = CM.Axis(f[1, 1], xlabel="Operator",
                 ylabel=L"|\langle \Omega|\, \mathcal{O}_i \,|n \rangle|^2",
                 xticks=(1:N_op, [L"\mathcal{O}_%$i" for i in 1:N_op]))

    # Plot histogram of overlaps
    color_grad = CM.cgrad(colors, n_levels)
    CM.barplot!(bar_xpos, bar_heights, stack=groups, color=color_grad[groups])

    # Add legend entry for each energy level
    elements = [CM.PolyElement(polycolor=color_grad[n]) for n in n_range]
    labels = [L"State $n=%$(i_n-1)$" for i_n in n_range]
    CM.Legend(f[1, 2], elements, labels)

    return ax
end
