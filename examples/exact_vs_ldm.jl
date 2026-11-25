using OpticalSpectrumModels
using CairoMakie

# `PointSample()` against `LDM()` on eight Voigt lines, one of them a pure
# Lorentzian (`sigma = 0`). LDM targets about 1% accuracy. It is not faster at
# this line count.

positions = [10.0, 22.0, 35.0, 48.0, 60.0, 71.0, 83.0, 92.0]
amplitudes = [1.0, 0.6, 1.4, 0.8, 1.2, 0.5, 1.0, 0.7]
sigmas = [0.30, 0.10, 0.45, 0.20, 0.15, 0.35, 0.0, 0.25]
gammas = [0.05, 0.25, 0.10, 0.15, 0.30, 0.05, 0.20, 0.10]

lines = [Line(positions[i], amplitudes[i], Voigt(sigmas[i], gammas[i])) for i in eachindex(positions)]

xs = range(0.0, 100.0; length=4000)

spectrum_exact = calculate_lines(xs, lines, PointSample())
spectrum_ldm = calculate_lines(xs, lines, LDM())

residual = spectrum_exact .- spectrum_ldm
max_rel_diff = maximum(abs, residual) / maximum(spectrum_exact)
l2_rel_diff = sqrt(sum(abs2, residual) / sum(abs2, spectrum_exact))
@info "PointSample vs LDM" max_rel_diff l2_rel_diff

fig = CairoMakie.Figure(size=(800, 500))

ax1 = CairoMakie.Axis(fig[1, 1]; ylabel="Intensity")
CairoMakie.lines!(ax1, xs, spectrum_exact; label="PointSample()", linewidth=2)
CairoMakie.lines!(ax1, xs, spectrum_ldm; label="LDM()", linewidth=2, linestyle=:dash)
CairoMakie.axislegend(ax1)

ax2 = CairoMakie.Axis(fig[2, 1]; xlabel="x", ylabel="PointSample − LDM")
CairoMakie.lines!(ax2, xs, residual; color=:red, linewidth=1)

pdf_path = joinpath(@__DIR__, "exact_vs_ldm.pdf")
CairoMakie.save(pdf_path, fig)
@info "PDF written" path=pdf_path
