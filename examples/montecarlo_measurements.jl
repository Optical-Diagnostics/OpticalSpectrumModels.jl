using SpectraUtils
using CairoMakie
using MonteCarloMeasurements
using Measurements

# Uncertain line parameters propagated through `calculate_lines`, with
# Measurements.jl (linear error propagation) and MonteCarloMeasurements.jl
# (particles). Each panel shows the mean ± one standard deviation. For the
# width the two differ: the peak height is nonlinear in it, and a few particles
# have a width near zero.

spectrum(xs, p) = calculate_lines(xs, [Line(0.0, p.amplitude, Gaussian(p.width))], PointSample()) .+ p.bg

xs = -5:0.1:5
nominal = (amplitude = 1.0, bg = 0.0, width = 1.0)

cases = [
    "Uncertain width"      => (:width, 1.0, 0.3),
    "Uncertain amplitude"  => (:amplitude, 1.0, 0.3),
    "Uncertain background" => (:bg, 0.0, 0.2),
]

band_measurements(ys) = (Measurements.value.(ys) .- Measurements.uncertainty.(ys),
                         Measurements.value.(ys) .+ Measurements.uncertainty.(ys))
band_particles(ys) = (pmean.(ys) .- pstd.(ys), pmean.(ys) .+ pstd.(ys))

fig = Figure(size = (1200, 450))
for (i, (title, (key, μ, σ))) in enumerate(cases)
    ax = Axis(fig[1, i]; title, xlabel = "x")
    band!(ax, xs, band_measurements(spectrum(xs, merge(nominal, (; key => Measurements.:±(μ, σ)))))...;
          alpha = 0.5, label = "Measurements")
    band!(ax, xs, band_particles(spectrum(xs, merge(nominal, (; key => MonteCarloMeasurements.:±(μ, σ)))))...;
          alpha = 0.5, label = "MonteCarloMeasurements")
    lines!(ax, xs, spectrum(xs, nominal); color = :black, label = "nominal")
    i == 1 && Legend(fig[2, 1:3], ax; orientation = :horizontal, framevisible = false)
end

pdf_path = joinpath(@__DIR__, "montecarlo_measurements.pdf")
save(pdf_path, fig)
@info "PDF written" path = pdf_path
