using OpticalSpectrumModels
using FittingUtils, NonlinearSolve, Turing
using ComponentArrays, Distributions, Random
using Statistics
using CairoMakie

# Fit one Lorentzian line to synthetic data, first with Levenberg-Marquardt,
# then with Turing (MCMC). Uses FittingUtils.jl.

Random.seed!(42)

# Synthetic data
spectrum(ws, p) = calculate_lines(ws, [Line(p.position, p.amplitude, Lorentzian(p.width))], PointSample())

true_params = ComponentVector(position=0.0, amplitude=1.0, width=0.5)
ws = -5:0.05:5
true_spec = spectrum(ws, true_params)

σ_noise = 0.1
noisy_spec = true_spec .+ rand(Normal(0, σ_noise), length(ws))

# FittingUtils takes the data as one `Normal` per point: observation and σ.
data = Normal.(noisy_spec, σ_noise)

# Levenberg-Marquardt fit
# The prior means are the initial guesses.
priors_lm = ComponentVector(position=Normal(0.1, 1.), amplitude=Normal(0.8, 0.5), width=Normal(0.5, 0.5))

println("Running LM Fit...")
@time res_lm = fit_model(spectrum, ws, data, priors_lm, LevenbergMarquardtFit())

println("\n--- Levenberg-Marquardt Results ---")
println("Coefficients: ", coef(res_lm))
println("Std Errors:   ", stderror(res_lm))

# 90% band of spectra drawn from the fitted parameter distributions. Each
# parameter is drawn from its own marginal, so correlations are ignored.
function spectrum_band(params; n = 500)
    samples = reduce(hcat, (spectrum(ws, rand.(params)) for _ in 1:n))
    lo = [quantile(row, 0.05) for row in eachrow(samples)]
    hi = [quantile(row, 0.95) for row in eachrow(samples)]
    return lo, hi
end

# Turing fit
priors_turing = ComponentVector(position=Normal(0.0, 1.0), amplitude=Normal(1.0, 0.5), width=Uniform(0.01, 2.0))

println("\nRunning Turing Fit...")
@time res_turing = fit_model(spectrum, ws, data, priors_turing, TuringFit(); num_samples=2000)

# `res_turing.marginals` holds the empirical posteriors, with the same
# coef/stderror interface as the LM result.
println("\n--- Turing Results ---")
println("Mean Params: ", coef(res_turing))

# Plot
fig = Figure(size=(1200, 500))
ax = CairoMakie.Axis(fig[1, 1], title="Fits with 90% bands", xlabel="Frequency", ylabel="Intensity")
scatter!(ax, ws, noisy_spec, label="Data", color=:gray, alpha=0.5)
lines!(ax, ws, true_spec, label="Truth", color=:black, linestyle=:dash)
for (res, name, color) in ((res_lm, "LM", :red), (res_turing, "Turing", :blue))
    band!(ax, ws, spectrum_band(res.marginals)..., color=(color, 0.2))
    lines!(ax, ws, spectrum(ws, coef(res)), label=name, color=color, linewidth=2)
end
axislegend(ax)

# Posteriors (Turing, filled) against the LM estimates (dashed).
gl = fig[1, 2] = GridLayout()
for (i, name) in enumerate(keys(res_turing.marginals))
    ax_p = CairoMakie.Axis(gl[i, 1], title=string(name), ylabel="Density")
    density!(ax_p, [rand(res_turing.marginals[name]) for _ in 1:2000], color=(:blue, 0.5), strokewidth=1)
    lines!(ax_p, res_lm.marginals[name], color=:red, linestyle=:dash)
    vlines!(ax_p, [true_params[name]], color=:black, linestyle=:dash)
end

pdf_path = joinpath(@__DIR__, "lorentzian_fit.pdf")
save(pdf_path, fig)
@info "PDF written" path=pdf_path
