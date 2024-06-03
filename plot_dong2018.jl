using Pkg
Pkg.add(["Revise", "Plots"])
using Revise
using LaTeXStrings
using Plots; gr()
(@isdefined Scattering_models) == false ? include("scattering_models.jl") : nothing
using .Scattering_models: scattering_dong2018


α = 10.
β_s = 0.03
# Impact parameter.
b = 2.4
# Observer plane coordinates to calculate the intensity.
u__array = vcat(LinRange(-6., 6., 300))
# Limits for root finding with lens equation.
u_min = -1.5*20
u_max = 1.5*20

I_result, I_err = scattering_dong2018(u__array, β_s, α, b, u_min, u_max)
title = L"Dong+2018 with $ α = %$α, β_s = %$β_s, b = %$b$ "

s = plot(u__array, I_result./I_result[1], xlabel="Epoch (normalized and centered)", ylabel="Normalized flux density",
         linewidth=3,
         size=(1000, 750), label="", title=title, xtickfontsize=16,
         ytickfontsize=16, xguidefontsize=18, yguidefontsize=18,
         titlefontsize=20, right_margin=5Plots.mm, left_margin=5Plots.mm,
         minorgrid=true, xlims=(-6., 6.), ylims=(0., 4.))
savefig(s, "Dong2018.png")
gui()
