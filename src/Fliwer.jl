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
export AbstractBody,Body,NoBody,sdf

#include("AutoBody.jl")
#export AutoBody,Bodies,measure,sdf,+,-

include("Capacity.jl")
export AbstractFluid,Fluid,measure!

include("Operators.jl")


include("Vizualize.jl")


end
