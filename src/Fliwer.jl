module Fliwer

using Reexport

@reexport using IterativeSolvers
@reexport using LinearAlgebra
@reexport using SparseArrays
@reexport using LsqFit
@reexport using Makie
@reexport using SpecialFunctions
@reexport using Test
@reexport using CairoMakie
@reexport using Roots
@reexport using DataFrames
@reexport using Polynomials
@reexport using StaticArrays
@reexport using ForwardDiff

include("Body.jl")
export AbstractBody,measure_sdf!

include("AutoBody.jl")
export AutoBody,Bodies,measure,sdf,+,-

include("Vizualize.jl")

include("Utils.jl")
export loc, δ, inside, @inside




end
