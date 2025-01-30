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
export CartesianSpaceTimeMesh

include("Body.jl")
export AbstractBody,Body,NoBody,sdf,+,⊖,c,measure, NoBody1D, NoBody2D, NoBody3D, curvature

include("Capacity.jl")
export AbstractCapacity,Capacity,measure!

include("Operators.jl")
export AbstractOperators, DiffusionOps, ConvectionOps, ẟ_m, δ_p, Σ_m, Σ_p, I
export AdvectionVecOps
export ∇ 
export SpaceTimeOps

include("Boundary.jl")
export AbstractBoundary, Dirichlet, Neumann, Robin, Periodic, AbstractInterfaceBC, ScalarJump, FluxJump, BorderConditions, InterfaceConditions

include("Phase.jl")
export Phase
export VectorPhase

include("Utils.jl")
export identify!, find_border, find_cut, eval_sdf 
export initialize_temperature_uniform!, initialize_temperature_square!, initialize_temperature_circle!, initialize_temperature_function!
export initialize_rotating_velocity_field, initialize_radial_velocity_field, initialize_poiseuille_velocity_field
export remove_zero_rows_cols!
export check_convergence

include("Solver.jl")
export TimeType, PhaseType, EquationType, Solver
export DiffusionSteadyMono, DiffusionSteadyDiph, DiffusionUnsteadyMono, DiffusionUnsteadyDiph
export solve_DiffusionSteadyMono!, solve_DiffusionSteadyDiph!, solve_DiffusionUnsteadyMono!, solve_DiffusionUnsteadyDiph!
export AdvectionUnsteadyMono
export solve_AdvectionUnsteadyMono!
export AdvectionDiffusionSteadyMono, AdvectionDiffusionSteadyDiph, AdvectionDiffusionUnsteadyMono, AdvectionDiffusionUnsteadyDiph
export solve_AdvectionDiffusionSteadyMono!, solve_AdvectionDiffusionSteadyDiph!, solve_AdvectionDiffusionUnsteadyMono!, solve_AdvectionDiffusionUnsteadyDiph!

include("SolverVec.jl")
export DiffusionVecUnsteadyMono, solve_DiffusionVecUnsteadyMono!
export ConvectionVecUnsteadyMono, solve_ConvectionVecUnsteadyMono!

include("MovingSolver.jl")
export MovingDiffusionUnsteadyMono, solve_MovingDiffusionUnsteadyMono!
export MovingDiffusionUnsteadyMono2, solve_MovingDiffusionUnsteadyMono2!
export MovingDiffusionUnsteadyDiph, solve_MovingDiffusionUnsteadyDiph!
export MovingDiffusionUnsteadyDiph2, solve_MovingDiffusionUnsteadyDiph2!

include("Solve.jl")
export solve!

include("Darcy.jl")
export DarcyFlow, solve_DarcyFlow!, solve_darcy_velocity
export DarcyFlowUnsteady, solve_DarcyFlowUnsteady!

include("Vizualize.jl")
export plot_solution, plot_profile, animate_solution
export plot_mesh, plot_body
export plot_solution_vector

include("Vtk.jl")
export write_vtk

end
