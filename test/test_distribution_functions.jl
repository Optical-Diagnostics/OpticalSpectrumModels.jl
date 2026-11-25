using OpticalSpectrumModels
using OpticalSpectrumModels: h, c0, kB
using Test
using ComponentArrays

const c2 = 100 * h * c0 / kB   # cm·K, matches Boltzmann's own exponent

@testset "DistributionFunction" begin
    s0 = State(:X, 0, nothing, 0.0, 1.0)
    s1 = State(:X, 1, nothing, 1000.0, 3.0)
    params = ComponentVector(T = 500.0, pop = 42.0)

    @testset "Boltzmann" begin
        Qfun(T) = 1.0 + 3.0 * exp(-c2 * 1000.0 / T)
        bd = Boltzmann(Qfun, 0.0, State(:X))
        @test bd(s1, params) ≈ 3.0 / Qfun(500.0) * exp(-1000.0 * c2 / 500.0)
        # Qfun sums over exactly these two states.
        @test bd(s0, params) + bd(s1, params) ≈ 1
        # A state outside `parent` contributes exactly zero.
        @test bd(State(:A, 0, nothing, 500.0, 1.0), params) == 0.0
        # Only a bare Boltzmann is LTE, not a wrapped one.
        @test is_LTE(bd)
        @test !is_LTE(ProductDistribution(bd))
    end

    @testset "StateResolved" begin
        dist = StateResolved(vibrational_state,
                             Dict(0 => ConstantDistribution(4.0), 1 => ConstantDistribution(6.0)))
        @test dist(s0, params) == 4.0
        @test dist(s1, params) == 6.0
        @test dist(State(:X, 2, nothing, 0.0, 1.0), params) == 0.0   # no entry: zero, not an error
        @test !is_LTE(dist)
    end

    @testset "ProductDistribution, ConstantDistribution, VariableDistribution" begin
        a, b, c = ConstantDistribution(2.0), ConstantDistribution(3.0), ConstantDistribution(5.0)
        @test ProductDistribution(a, b, c)(s0, params) == 30.0
        @test ProductDistribution((a, b))(s0, params) == 6.0
        @test a(s0) == 2.0   # params is optional
        @test VariableDistribution(:pop)(s0, params) == 42.0

        @test DistributionFunction(2.5) isa ConstantDistribution
        vd = VariableDistribution(:x)
        @test DistributionFunction(vd) === vd
    end
end
