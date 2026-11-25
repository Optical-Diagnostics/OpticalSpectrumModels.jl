using OpticalSpectrumModels
using Test
using ComponentArrays
using ForwardDiff

include("synthetic.jl")

const NL = (number_density_function = (p) -> p[:N],
            path_length_function    = (p) -> p[:L])

@testset "RadianceSpectrum" begin

    @testset "LTE fast path ≡ general path (Kirchhoff identity)" begin
        # For one transition, j/α is exactly planck_radiance(ν̃_t, T) at the rest
        # wavenumber. The fast path uses planck_radiance(ν, T) at the grid point.
        # The two agree when lines are narrow, as here.
        mk(df) = RadianceSpectrum(;
            species = SpeciesModel(; transitions, distribution_function = df,
                                     lineshape_function = (t, p) -> Lorentzian(0.01),
                                     number_density_function = (p) -> p[:N],
                                     method = PointSample(wing_cut = 20.0)),   # ±0.2 cm⁻¹
            path_length_function = (p) -> p[:L])
        xs = sort(vcat(range(999.7, 1000.3; length = 61), range(1049.7, 1050.3; length = 61)))
        params = ComponentVector(T = 800.0, N = 1e15, L = 10.0)
        @test mk(bd)(xs, params) ≈ mk(ProductDistribution(bd))(xs, params) rtol = 1e-3
    end

    @testset "optically-thin limit ≈ manual j·L" begin
        # The general (α, j) path, forced by wrapping the Boltzmann, reduces
        # to j·L at small αL.
        bd10 = Boltzmann(T -> 10.0, 0.0, State(:X))
        kw = (transitions, distribution_function = ProductDistribution(bd10),
              lineshape_function = (t, p) -> Lorentzian(0.5), NL...)
        xs = collect(range(950.0, 1100.0; length = 300))
        N, L, T = 2.5, 3.0, 800.0
        params = ComponentVector(T = T, N = N, L = L)

        manual = zeros(length(xs))
        for t in transitions
            f_u = t.upper.degeneracy / 10.0 * exp(-c2 * t.upper.energy / T)
            manual .+= N * L * (h * c0_cm * t.wavenumber / 4π) * t.A * f_u .*
                       Lorentzian(0.5).(xs .- t.wavenumber)
        end
        @test RadianceSpectrum(; kw...)(xs, params) ≈ manual rtol = 1e-6

        # EmissionSpectrum with amplitude = SpectralRadiance() is exactly this
        # limit; the PhotonRate default differs per line by h·c·ν̃/4π.
        @test EmissionSpectrum(; kw..., amplitude = SpectralRadiance())(xs, params) ≈ manual rtol = 1e-6
        em1_photon = EmissionSpectrum(; kw..., transitions = [t1])
        em1_energy = EmissionSpectrum(; kw..., transitions = [t1], amplitude = SpectralRadiance())
        @test em1_energy(xs, params) ≈ em1_photon(xs, params) .* (h * c0_cm * t1.wavenumber / 4π)
    end

    @testset "non-equilibrium (StateResolved per-level population)" begin
        # Each vibrational level has its own partition sum.
        Qv = Dict(0 => 4.0, 1 => 6.0)
        df = StateResolved(vibrational_state, Dict(
            0 => Boltzmann(T -> 4.0, 0.0, State(:X, 0)),
            1 => Boltzmann(T -> 6.0, 0.0, State(:X, 1))))
        sp = RadianceSpectrum(; transitions, distribution_function = df,
                                lineshape_function = (t, p) -> Lorentzian(0.5), NL...)
        xs = collect(range(950.0, 1100.0; length = 300))
        N, L, T = 1e15, 10.0, 800.0
        params = ComponentVector(T = T, N = N, L = L)

        # α and j by hand, then the slab solution I = j/α·(1 − e^{−αL}).
        f(s) = s.degeneracy / Qv[s.vibrational_state] * exp(-c2 * s.energy / T)
        α = zeros(length(xs))
        j = zeros(length(xs))
        for t in transitions
            φ  = Lorentzian(0.5).(xs .- t.wavenumber)
            gl, gu = t.lower.degeneracy, t.upper.degeneracy
            σ0 = gu / gl * t.A / (8π * c0_cm * t.wavenumber^2)
            α .+= N * σ0 * (f(t.lower) - f(t.upper) * gl / gu) .* φ
            j .+= N * (h * c0_cm * t.wavenumber / 4π) * t.A * f(t.upper) .* φ
        end
        @test sp(xs, params) ≈ @.(j / α * -expm1(-α * L)) rtol = 1e-10
    end

    @testset "BinAverage converges to the true bin-averaged radiance" begin
        # The RTE is nonlinear in α, so combining bin-averaged α and j is
        # second-order accurate in the bin width.
        bd10 = Boltzmann(T -> 10.0, 0.0, State(:X))
        mk(m) = RadianceSpectrum(;
            species = SpeciesModel(; transitions = [t1],
                                     distribution_function = ProductDistribution(bd10),  # general α/j path
                                     lineshape_function = (t, p) -> Lorentzian(0.05),
                                     number_density_function = (p) -> p[:N], method = m),
            path_length_function = (p) -> p[:L])
        params = ComponentVector(T = 1500.0, N = 1.0e8, L = 100.0)   # τ_peak ≈ 1.6

        # The reference grid must be much finer than the finest bin, or its
        # own error sets the floor.
        xf = collect(range(999.0, 1001.0; length = 200001))
        yf = mk(PointSample())(xf, params)
        errs = map((0.2, 0.05, 0.01)) do Δ
            xc = collect(range(999.0, 1001.0; step = Δ))
            yb = mk(BinAverage(bin_width = Δ))(xc, params)
            ref = map(xc) do x
                i0 = searchsortedfirst(xf, x - Δ/2)
                i1 = searchsortedlast(xf, x + Δ/2)
                sum(@view yf[i0:i1]) / (i1 - i0 + 1)
            end
            maximum(abs, yb .- ref) / maximum(ref)
        end
        @test issorted(errs; rev = true)   # strictly converging
        @test errs[end] < 5.0e-3           # and small once the bin is sub-linewidth
    end

    @testset "AD" begin
        sp = RadianceSpectrum(;
            species = SpeciesModel(; transitions, distribution_function = bd,
                                     lineshape_function = (t, p) -> Lorentzian(0.5),
                                     number_density_function = (p) -> p[:N],
                                     method = PointSample(wing_cut = 100.0)),
            path_length_function = (p) -> p[:L])
        xs = collect(range(950.0, 1100.0; length = 50))
        f(p) = sp(xs, p)
        J = ForwardDiff.jacobian(f, ComponentVector(T = 800.0, N = 1e15, L = 10.0))

        # τ depends on N and L only through N·L, so N·∂/∂N = L·∂/∂L.
        @test 1e15 .* J[:, 2] ≈ 10.0 .* J[:, 3] rtol = 1e-10
        dT = 1e-3
        fd = (f(ComponentVector(T = 800.0 + dT, N = 1e15, L = 10.0)) .-
              f(ComponentVector(T = 800.0 - dT, N = 1e15, L = 10.0))) ./ (2dT)
        @test J[:, 1] ≈ fd rtol = 1e-6
    end

    @testset "background (instrumental) vs incident_radiance (physical)" begin
        mk(df; kw...) = RadianceSpectrum(;
            species = SpeciesModel(; transitions, distribution_function = df,
                                     lineshape_function = (t, p) -> Lorentzian(0.05),
                                     number_density_function = (p) -> p[:N],
                                     method = PointSample(wing_cut = 100.0)),   # ±5 cm⁻¹ per line
            path_length_function = (p) -> p[:L], kw...)
        xs     = collect(range(950.0, 1200.0; length = 200))   # 1200: α = 0 exactly
        params = ComponentVector(T = 800.0, N = 1e19, L = 100.0)  # optically thick at line centres

        r0 = mk(bd)(xs, params)
        # Exactly zero, not NaN, where no line contributes and I₀ = 0.
        @test all(isfinite, r0)
        @test r0[end] == 0.0

        # Instrumental background: added everywhere, never attenuated by the
        # slab, on both the LTE and the general path.
        r_bg = mk(bd; background = 0.7)(xs, params)
        @test r_bg ≈ r0 .+ 0.7
        pd = ProductDistribution(bd)
        @test mk(pd; background = 0.7)(xs, params) ≈ mk(pd)(xs, params) .+ 0.7
        ret = similar(r_bg); mk(bd; background = 0.7)(ret, xs, params)
        @test ret ≈ r_bg

        # Incident radiance: fully transmitted where α = 0, absorbed at an
        # optically-thick line centre.
        r_inc = mk(bd; incident_radiance = 2.0)(xs, params)
        @test r_inc[end] ≈ r0[end] + 2.0
        i_c = argmin(abs.(xs .- 1000.0))
        @test r_inc[i_c] - r0[i_c] < 1e-6
    end

    @testset "_combine_radiance: α → 0 values and derivatives" begin
        cr = OpticalSpectrumModels._combine_radiance
        j, L, I₀ = 2.0, 5.0, 0.3
        # Exact α = 0: optically-thin limit j·L plus untouched I₀.
        @test cr(j, 0.0, L, I₀) == j * L + I₀
        # Accuracy across the series/expm1 switch (|τ| < 1e-4), BigFloat reference.
        @test all((1e-8, 1e-5, 9.9e-5, 1.01e-4, 1e-3, 1.0, 10.0)) do τ
            ref = j * L * Float64((1 - exp(-big(τ))) / big(τ)) + I₀ * exp(-τ)
            isapprox(cr(j, τ / L, L, I₀), ref; rtol = 1e-13)
        end
        # d/dα at α = 0 is −j·L²/2 − I₀·L.
        @test ForwardDiff.derivative(α -> cr(j, α, L, I₀), 0.0) ≈ -j * L^2 / 2 - I₀ * L
    end

    @testset "in-place form: parity, buffer reuse, allocations" begin
        xs     = collect(range(950.0, 1100.0; length = 300))
        params = ComponentVector(T = 800.0, N = 1e15, L = 10.0)
        for df in (bd, ProductDistribution(bd))   # LTE fast path and general path
            sp = RadianceSpectrum(; transitions, distribution_function = df,
                                    lineshape_function = (t, p) -> Lorentzian(0.5), NL...)
            ref = sp(xs, params)
            ret = similar(ref)
            sp(ret, xs, params)
            sp(ret, xs, params)   # reuse must not corrupt
            @test ret == ref

            # In-place allocates less than allocating, and no more on a 10x grid.
            a_inplace = @allocated sp(ret, xs, params)
            @test a_inplace < @allocated sp(xs, params)
            xs_big  = collect(range(950.0, 1100.0; length = 3000))
            ret_big = sp(xs_big, params)
            sp(ret_big, xs_big, params)
            @test (@allocated sp(ret_big, xs_big, params)) <= a_inplace + 512
        end
    end

    @testset "rejected inputs" begin
        xs     = collect(range(950.0, 1100.0; length = 50))
        params = ComponentVector(T = 800.0, N = 1e15, L = 10.0)

        # Convolving α with the instrument before 1 − e^{−αL} is wrong at
        # τ ≳ 1, so both paths refuse a per-line instrument.
        inst = (transitions, lineshape_function = (t, p) -> Lorentzian(0.05),
                instrument_function = (wn, p) -> Gaussian(0.02), NL...)
        @test_throws ArgumentError RadianceSpectrum(; inst..., distribution_function = bd)(xs, params)
        @test_throws ArgumentError RadianceSpectrum(; inst..., distribution_function = ProductDistribution(bd))(xs, params)

        # A Boltzmann parent that excludes the lower states is rejected. One
        # that excludes the upper states is allowed.
        tA = [ElectricDipoleTransition(State(:X, nothing, nothing, 0.0, 1.0),
                                       State(:A, nothing, nothing, 30000.0, 3.0), 30000.0, 1.0e6)]
        ax = (transitions = tA, lineshape_function = (t, p) -> Lorentzian(0.05), NL...)
        @test_throws ArgumentError RadianceSpectrum(; ax...,
            distribution_function = Boltzmann(T -> 3.0, 30000.0, State(:A)))
        sp_x = RadianceSpectrum(; ax..., distribution_function = Boltzmann(T -> 1.0, 0.0, State(:X)))
        @test all(>=(0), sp_x(collect(29999.0:0.5:30001.0), params))

        empty_tr = ElectricDipoleTransition[]
        @test_throws ArgumentError EmissionSpectrum(;
            transitions = empty_tr, distribution_function = bd,
            lineshape_function = (t, p) -> Lorentzian(0.05))(xs, params)
        @test_throws ArgumentError RadianceSpectrum(;
            transitions = empty_tr, distribution_function = bd,
            lineshape_function = (t, p) -> Lorentzian(0.05), NL...)(xs, params)
        @test_throws ArgumentError calculate_radiance_lines(xs, OpticalSpectrumModels._DualLine[], PointSample())

        for method in (PointSample(min_amplitude = 1.0e-30), BinAverage(bin_width = 0.1, min_amplitude = 1.0e-30))
            sp = RadianceSpectrum(; species = SpeciesModel(; transitions = [t1],
                    distribution_function = ConstantDistribution(0.5),
                    lineshape_function = (t, p) -> Lorentzian(0.1), method),
                path_length_function = p -> 1.0)
            @test_throws ArgumentError sp(xs, params)
        end
    end

    @testset "non-LTE path differentiates with every shape" begin
        # Dual amplitudes over Float64 widths: the shape is converted to the
        # dual type, also where the numeric parameter is not the first one.
        Qfun(T) = 1.0 + 3.0 * exp(-c2 * 1000.0 / T) + 5.0 * exp(-c2 * 1050.0 / T)
        nonlte = ProductDistribution(Boltzmann(Qfun, 0.0, State(:X)))
        xs = range(990.0, 1060.0, length = 701)
        for shape in (Voigt(0.05, 0.02), PseudoVoigt(0.05, 0.02), Gaussian(0.05), Triangular(0.1))
            sp = RadianceSpectrum(; transitions, distribution_function = nonlte,
                                    lineshape_function = (t, p) -> shape,
                                    number_density_function = p -> 1e16, path_length_function = p -> 1.0)
            f(T) = sum(sp(xs, ComponentVector(T = T)))
            @test ForwardDiff.derivative(f, 1000.0) ≈ (f(1000.5) - f(999.5)) rtol = 1e-5
        end
    end

    @testset "calculate_radiance_lines! checks its buffers and grid" begin
        Qfun(T) = 1.0 + 3.0 * exp(-c2 * 1000.0 / T) + 5.0 * exp(-c2 * 1050.0 / T)
        bd = Boltzmann(Qfun, 0.0, State(:X))
        p = ComponentVector(T = 1000.0)
        dls = [radiance_line(t, bd, (t, p) -> Lorentzian(0.1), p) for t in transitions]
        xs = range(990.0, 1060.0, length = 71)
        for method in (PointSample(), BinAverage(bin_width = 1.0))
            @test_throws DimensionMismatch calculate_radiance_lines!(zeros(3), zeros(71), xs, dls, method)
            @test_throws DimensionMismatch calculate_radiance_lines!(zeros(71), zeros(3), xs, dls, method)
        end
        for method in (PointSample(wing_cut = 10.0), BinAverage(bin_width = 1.0, wing_cut = 10.0))
            @test_throws ArgumentError calculate_radiance_lines(reverse(xs), dls, method)
        end
    end
end
