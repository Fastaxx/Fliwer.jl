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
export CartesianMesh, nC

include("Body.jl")
export AbstractBody,Body,NoBody,sdf,+,⊖,c,measure

include("Capacity.jl")
export AbstractCapacity,Capacity,measure!

include("Operators.jl")
export AbstractOperators, DiffusionOps, ConvectionOps, ẟ_m, δ_p, Σ_m, Σ_p, I

include("Boundary.jl")
export AbstractBoundary, Dirichlet, Neumann, Robin, Periodic, AbstractInterface, ScalarJump, FluxJump

include("Utils.jl")
export find_border, eval_sdf, find_cut

include("Solver.jl")
export TimeType, PhaseType, EquationType, Solver, DiffusionSolver, AdvectionSolver, AdvectionDiffusionSolver, initialize!, solve!


end
