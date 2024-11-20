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
using GLMakie
using StaticArrays
using ForwardDiff
using CartesianGeometry
using WriteVTK

include("Mesh.jl")
export CartesianMesh, nC, MeshTag

include("Body.jl")
export AbstractBody,Body,NoBody,sdf,+,⊖,c,measure

include("Capacity.jl")
export AbstractCapacity,Capacity,measure!

include("Operators.jl")
export AbstractOperators, DiffusionOps, ConvectionOps, ẟ_m, δ_p, Σ_m, Σ_p, I

include("Boundary.jl")
export AbstractBoundary, Dirichlet, Neumann, Robin, Periodic, AbstractInterfaceBC, ScalarJump, FluxJump, BorderConditions, InterfaceConditions

include("Phase.jl")
export Phase

include("Utils.jl")
export identify!, find_border, find_cut, eval_sdf

include("Solver.jl")
export TimeType, PhaseType, EquationType, Solver, DiffusionSteadyMono, DiffusionSteadyDiph, DiffusionUnsteadyMono, DiffusionUnsteadyDiph, AdvectionSteadyMono, solve!

include("Vizualize.jl")
export plot_solution, plot_profile, animate_solution

include("Vtk.jl")
export write_vtk

end
