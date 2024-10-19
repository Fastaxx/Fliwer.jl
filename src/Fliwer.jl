module Fliwer

using Reexport

using IterativeSolvers
using LinearAlgebra
using SparseArrays
using LsqFit
using Makie
using SpecialFunctions
using Test
using CairoMakie
using Roots
using DataFrames
using Polynomials
using StaticArrays
using ForwardDiff

include("Capacity.jl")
export AbstractFluid,Fluid,evaluate_capacity

include("Body.jl")
export AbstractBody,measure_sdf!

include("AutoBody.jl")
export AutoBody,Bodies,measure,sdf,+,-

include("Vizualize.jl")


end
