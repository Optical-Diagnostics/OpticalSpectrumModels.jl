using Test, JET
using OpticalSpectrumModels
using ComponentArrays

const TM = (OpticalSpectrumModels,)   # target_modules shorthand

@testset "JET call analysis" begin
    params = ComponentVector(T = 500.0)
    s  = State(:X, 0, 1)
    vd = VariableDistribution(:T)
    @test_call target_modules=TM ProductDistribution(vd, vd)(s, params)
    @test_call target_modules=TM StateResolved(electronic_state, Dict(:X => vd))(s, params)

    transitions = [ElectricDipoleTransition(State(:X, 0, i), State(:A, 1, i), 1000.0 + float(i), 1e6)
                   for i in 1:5]
    sp = EmissionSpectrum(; transitions, distribution_function = vd,
                            lineshape_function = (t, p) -> Lorentzian(0.1))
    xs = collect(range(999.0, 1007.0, length = 50))
    @test_call target_modules=TM sp(xs, params)
    @test_call target_modules=TM sp(zeros(length(xs)), xs, params)
end
