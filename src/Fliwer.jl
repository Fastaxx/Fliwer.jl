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
using StaticArrays
using ForwardDiff
using CartesianGeometry

include("Mesh.jl")
export CartesianMesh, nodes, centers, nC

include("Body.jl")
export AbstractBody,Body,NoBody,sdf,+,⊖,c,measure

include("Capacity.jl")
export AbstractFluid,Fluid,measure!

include("Operators.jl")
export  AbstractOperators, DiffusionOps, ConvectionOps, ẟ_m, δ_p, Σ_m, Σ_p, I

include("Vizualize.jl")


end
