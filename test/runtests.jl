using SafeTestsets

@safetestset "convolve_spectrum" begin
    include("test_convolve_spectrum.jl")
end

@safetestset "Linewidths and broadening" begin
    include("test_linewidths.jl")
end

@safetestset "ChebyshevBackground" begin
    include("test_backgrounds.jl")
end

@safetestset "EmissionSpectrum" begin
    include("test_emission.jl")
end

@safetestset "HITRAN submodule" begin
    include("test_hitran.jl")
end

@safetestset "LineCache / fill_lines!" begin
    include("test_line_cache.jl")
end

@safetestset "RadianceSpectrum" begin
    include("test_radiance.jl")
end

@safetestset "Absorption-side models" begin
    include("test_absorption_models.jl")
end

@safetestset "Multi-species models" begin
    include("test_multispecies.jl")
end

@safetestset "Multi-slab models" begin
    include("test_multislab.jl")
end

@safetestset "StaticLineshapes" begin
    include("test_static_lineshapes.jl")
end

@safetestset "DistributionFunction" begin
    include("test_distribution_functions.jl")
end

@safetestset "ExoMol extension" begin
    include("test_exomol.jl")
end

@safetestset "Code quality (Aqua.jl)" begin
    using Aqua, OpticalSpectrumModels
    # DataInterpolations is used only by the ExoMol extension, which can load
    # only the parent package's deps.
    Aqua.test_all(OpticalSpectrumModels; stale_deps = (ignore = [:DataInterpolations],))
end

@safetestset "JET analysis" begin
    include("test_jet.jl")
end
