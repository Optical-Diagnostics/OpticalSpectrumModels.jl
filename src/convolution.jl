# SpectralLines' own grid checks and full linear convolution: FFT for hardware
# floats, a direct loop for anything else (ForwardDiff duals).
import SpectralLines: check_ascending, uniform_step, conv_full

"""
    convolve_spectrum(model, xs, params, instrument_function) -> Vector

FFT-convolve any `model(xs, params) -> Vector` with `instrument_function`.
`model` is evaluated on a grid padded by
`convolution_support(instrument_function)` on each side, so the crop back to `xs`
avoids edge effects (exactly `±fwhm` for `Triangular`, zero beyond that).

Needed for instrument shapes with no closed-form combination with the line shape,
e.g. `Triangular`; `Gaussian`/`Lorentzian`/`Voigt` combine analytically and are
applied per line inside the models. Pass `instrument_function = (wn, p) ->
nothing` to `model` so its per-line convolution is a no-op, and make `xs` fine
enough to resolve the natural linewidth.

`xs` must be ascending and uniformly spaced; anything else throws an
`ArgumentError`. Spacings within 1% of the first one count as uniform.
"""
function convolve_spectrum(model, xs, params, instrument_function)
    check_ascending(xs, "convolve_spectrum")
    dx = uniform_step(xs)
    lo, hi = convolution_support(instrument_function)
    n_pad = ceil(Int, max(-lo, hi) / dx)
    xs_padded = _padded_grid(xs, n_pad, dx)
    spectrum_padded = model(xs_padded, params)

    n_kernel = 2n_pad + 1
    kernel = instrument_function.(range(-n_pad * dx, n_pad * dx, length=n_kernel)) .* dx
    kernel ./= sum(kernel)

    return conv_full(spectrum_padded, kernel)[2n_pad+1:end-2n_pad]
end

# `xs` extended by `n_pad` points of spacing `dx` on each side. A range stays a
# range, which `LDM` needs. Any other grid keeps its own points between the
# pads.
_padded_grid(xs::AbstractRange, n_pad, dx) =
    range(first(xs) - n_pad * dx; step=dx, length=length(xs) + 2n_pad)
_padded_grid(xs, n_pad, dx) =
    vcat(range(first(xs) - n_pad * dx; step=dx, length=n_pad), xs,
         range(last(xs) + dx; step=dx, length=n_pad))
