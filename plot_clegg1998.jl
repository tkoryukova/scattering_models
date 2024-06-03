using Pkg
Pkg.add(["Revise", "Plots"])
using Revise
using LaTeXStrings
using Plots; gr()
(@isdefined Scattering_models) == false ? include("scattering_models.jl") : nothing
using .Scattering_models: scattering_clegg1998


α = 25.
β_s = 1.0
# Observer plane coordinates to calculate the intensity.
u__array = vcat(LinRange(-20., 20., 300))
# Limits for root finding with lens equation. Must be at least as wide as range of ``u_``.
u_min = -30.
u_max = 30.

I_result, I_err = scattering_clegg1998(u__array, β_s, α, u_min, u_max)
title = L"Clegg+1998 with $ α = %$α, β_s = %$β_s$"
s = plot(u__array, I_result, xlabel="Epoch (normalized and centered)", ylabel="Normalized flux density",
         linewidth=3, size=(1000, 750), label="", title=title, xtickfontsize=16,
         ytickfontsize=16, xguidefontsize=18, yguidefontsize=18,
         titlefontsize=20, right_margin=5Plots.mm, left_margin=5Plots.mm,
         minorgrid=true, xlims=(-20., 20.), ylims=(0., 2.))
savefig(s, "Clegg1998.png")
gui()
