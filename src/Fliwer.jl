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

include("Mesh.jl")
export CartesianMesh, nodes, centers, nC

include("Body.jl")
export AbstractBody,measure_sdf!

include("AutoBody.jl")
export AutoBody,Bodies,measure,sdf,+,-

include("Capacity.jl")
export AbstractFluid,Fluid,evaluate_capacity

include("Vizualize.jl")


end
