using QuantumTimeSteppers
using Test

@testset "QuantumTimeSteppers.jl" begin
    H = HermiteFSE()
    f = x -> (1+x^2)*exp(-x^2 / 2)
    F = Fun(f,H)
    x = 0.78
    @test F(x) ≈ f(x)

    f = x -> cos(10x)*exp(-x^2)
    F = Fun(f,H)
    x = 0.78
    @test F(x) ≈ f(x)

    v = [0.83;0.78;0.65;1.11]
    trans = QuantumTimeSteppers.plan_transform(H,v)
    itrans = QuantumTimeSteppers.plan_itransform(H,v)
    @test v ≈ itrans*(trans*v)
    @test v ≈ trans*(itrans*v)


    H = HermiteFSE()
    f = x -> sech(x)^2
    fprime = x -> -2tanh(x)*sech(x)^2
    F = Fun(f,H)
    Fprime = differentiate(F)
    x = 0.78
    @test Fprime(x) ≈ fprime(x)
    @test (Derivative(H,2)*F)(x) ≈ (differentiate(Fprime))(x)

end
