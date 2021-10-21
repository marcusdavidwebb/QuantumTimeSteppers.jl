module QuantumTimeSteppers

using ApproxFun, FastTransforms, Plots, LinearAlgebra, FastGaussQuadrature
import ApproxFun: points, plan_transform, plan_itransform, TransformPlan, ITransformPlan, domain, canonicaldomain, spacescompatible, evaluate
import ApproxFun: Multiplication, ConcreteMultiplication, Derivative, ConcreteDerivative, space, rangespace, bandwidths
import Base: *, first, last, getindex
export HermiteFSE
export Lie_Trotter_evolve, Strang_evolve, Zassenhaus2_evolve

include("splitting_methods.jl")
include("HermiteFSE.jl")

end
