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
export initialize_temperature_uniform!, initialize_temperature_square!, initialize_temperature_circle!, initialize_temperature_function!
export initialize_rotating_velocity_field, initialize_radial_velocity_field, initialize_poiseuille_velocity_field

include("Solver.jl")
export TimeType, PhaseType, EquationType, Solver
export DiffusionSteadyMono, DiffusionSteadyDiph, DiffusionUnsteadyMono, DiffusionUnsteadyDiph
export solve_DiffusionSteadyMono!, solve_DiffusionSteadyDiph!, solve_DiffusionUnsteadyMono!, solve_DiffusionUnsteadyDiph!
export AdvectionUnsteadyMono
export solve_AdvectionUnsteadyMono!
export AdvectionDiffusionSteadyMono, AdvectionDiffusionSteadyDiph, AdvectionDiffusionUnsteadyMono, AdvectionDiffusionUnsteadyDiph
export solve_AdvectionDiffusionSteadyMono!, solve_AdvectionDiffusionSteadyDiph!, solve_AdvectionDiffusionUnsteadyMono!, solve_AdvectionDiffusionUnsteadyDiph!

include("Solve.jl")
export solve!

include("Vizualize.jl")
export plot_solution, plot_profile, animate_solution

include("Vtk.jl")
export write_vtk

end
