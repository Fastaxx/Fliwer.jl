mutable struct VectorSolver{TT<:TimeType, PT<:PhaseType, ET<:EquationType}
    time_type::TT
    phase_type::PT
    equation_type::ET
    A::Union{SparseMatrixCSC{Float64, Int}, Nothing}
    b::Union{Vector{Float64}, Nothing}
    x::Union{Vector{Float64}, Nothing}
    ch::IterativeSolvers.ConvergenceHistory
    states::Vector{Any}
end

# Navier-Stokes Monophasic Unsteady
function NavierStokesUnsteadyMono(velocity, bc_b, Δt, Tend, x0)
    println("Création du solveur Navier-Stokes Monophasic Unsteady")
    println("- Monophasic problem")
    println("- Unsteady problem")

    s = VectorSolver(Unsteady, Monophasic, DiffusionAdvection, nothing, nothing, nothing, ConvergenceHistory(), [])

    if typeof(velocity)==Velocity{1}
        s.A = build_navier_stokes_matrix_init(velocity.operator, velocity.capacities, velocity.ρ, velocity.Re, Δt, x0)
        s.b = build_navier_stokes_rhs_init(velocity.operator, velocity.capacities, velocity.source, velocity.ρ, velocity.Re, Δt, x0)
        push!(s.states, x0)
    elseif typeof(velocity)==Velocity{2}
        println("2D")
    end

    # BC_border_mono!(s.A, s.b, bc_b, phase.capacity.mesh)
    return s
end

function build_navier_stokes_matrix_init(operator, capacity_u, ρ, Re, Δt, x0)
    n = prod(operator.size)

    #A = spzeros(Float64, n, n)
    Vx = operator.V[1]
    Gx = operator.G[1]
    Wꜝx= operator.Wꜝ[1]
    Hx = operator.H[1]
    Cx = operator.C[1]
    Kx = operator.K[1]

    # Compute the blocks
    block1 = Vx/Δt + 1/(2*Re) * Gx' * Wꜝx * Gx + 1/2 * Cx + 1/4 * Kx # Accumulation + Diffusion + Convection (t=0)
    block2 = 1/(2*Re) * Gx' * Wꜝx * Hx + 1/4 * Kx # Diffusion + Convection (t=0)
    block3 = 1/ρ * Wꜝx * Gx
    block4 = 1/ρ * Wꜝx * Hx
    block5 = I(n)
    block6 = -(Gx' + Hx')
    block7 = Hx'

    # Fill the matrix
    A = [block1 block2 block3 block4;
         zeros(n, n) block5 zeros(n, n) zeros(n, n);
         block6 block7 zeros(n, n) zeros(n, n);
         zeros(n, n) zeros(n, n) zeros(n, n) zeros(n, n)]

    return A
end

function build_navier_stokes_rhs_init(operator, capacity_u, source, ρ, Re, Δt, x0)
    n = prod(operator.size)
    b = zeros(4n)

    Vx = operator.V[1]
    Gx = operator.G[1]
    Wꜝx= operator.Wꜝ[1]
    Hx = operator.H[1]
    Cx = operator.C[1]
    Kx = operator.K[1]

    uₒ = (x0[1:n],)
    uᵧ = x0[n+1:2n]

    fₒn, fₒn1 = build_source(operator, source, 0, capacity_u[1]), build_source(operator, source, Δt, capacity_u[1])

    # Compute the convective terms
    convx = -1/2 * Cx * uₒ[1] - 1/4 * Kx * uₒ[1] - 1/2 * Kx * uᵧ

    b1 = Vx/Δt * uₒ[1] - 1/(2*Re) * Gx' * Wꜝx * Gx * uₒ[1] - 1/(2*Re) * Gx' * Wꜝx * Hx * uᵧ + convx + Vx/2 * (fₒn + fₒn1) # Accumulation + Diffusion + Convection (t=0) + Source
    b2 = zeros(n)
    b3 = zeros(n)
    b4 = zeros(n)

    b = [b1; b2; b3; b4]
    return b
end

function build_navier_stokes_matrix(operator, capacity_u, ρ, Re, Δt, xn, xn1)
    # xn refer to xn and xn1 refer to xn-1
    n = prod(operator.size)

    uₒn = (xn[1:n],)
    uᵧn = xn[n+1:2n]
    uₒn1 = (xn1[1:n],)
    uᵧn1 = xn1[n+1:2n]

    Vx = operator.V[1]
    Gx = operator.G[1]
    Wꜝx= operator.Wꜝ[1]
    Hx = operator.H[1]
    Cxn = δ_p(n) * spdiagm(0 => (Σ_m(n) * capacity_u[1].A[1] * uₒn[1])) * Σ_m(n)
    Cxn1 = δ_p(n) * spdiagm(0 => (Σ_m(n) * capacity_u[1].A[1] * uₒn1[1])) * Σ_m(n)
    Kxn = spdiagm(0 => Σ_p(n) * Hx' * uᵧn)
    Kxn1 = spdiagm(0 => Σ_p(n) * Hx' * uᵧn1)

    # Compute the blocks
    block1 = Vx/Δt + 1/(2*Re) * Gx' * Wꜝx * Gx + 3/4 * Cxn + 3/8 * Kxn  - 1/4 * Cxn1 - 1/8 * Kxn1 # Accumulation + Diffusion + Convection
    block2 = 1/(2*Re) * Gx' * Wꜝx * Hx + 3/8 * Kxn - 1/8 * Kxn1 # Diffusion + Convection
    block3 = 1/ρ * Wꜝx * Gx
    block4 = 1/ρ * Wꜝx * Hx
    block5 = I(n)
    block6 = -(Gx' + Hx')
    block7 = Hx'

    # Fill the matrix
    A = [block1 block2 block3 block4;
         zeros(n, n) block5 zeros(n, n) zeros(n, n);
         block6 block7 zeros(n, n) zeros(n, n);
         zeros(n, n) zeros(n, n) zeros(n, n) zeros(n, n)]

    return A
end

function build_navier_stokes_rhs(operator, capacity_u, source, t, ρ, Re, Δt, xn, xn1)
    n = prod(operator.size)
    b = zeros(4n)

    uₒn = (xn[1:n],)
    uᵧn = xn[n+1:2n]
    uₒn1 = (xn1[1:n],)
    uᵧn1 = xn1[n+1:2n]

    fₒn, fₒn1 = build_source(operator, source, t, capacity_u[1]), build_source(operator, source, t+Δt, capacity_u[1])

    Vx = operator.V[1]
    Gx = operator.G[1]
    Wꜝx= operator.Wꜝ[1]
    Hx = operator.H[1]
    Cxn = δ_p(n) * spdiagm(0 => (Σ_m(n) * capacity_u[1].A[1] * uₒn[1])) * Σ_m(n)
    Cxn1 = δ_p(n) * spdiagm(0 => (Σ_m(n) * capacity_u[1].A[1] * uₒn1[1])) * Σ_m(n)
    Kxn = spdiagm(0 => Σ_p(n) * Hx' * uᵧn)
    Kxn1 = spdiagm(0 => Σ_p(n) * Hx' * uᵧn1)

    # Compute the convective terms
    convxn = -3/4 * Cxn * uₒn[1] - 3/8 * Kxn * uₒn[1] - 3/8 * Kxn * uᵧn
    convxn1 = 1/4 * Cxn1 * uₒn1[1] + 1/8 * Kxn1 * uₒn1[1] + 1/8 * Kxn1 * uᵧn1

    b1 = Vx/Δt * uₒn[1] - 1/(2*Re) * Gx' * Wꜝx * Gx * uₒn[1] - 1/(2*Re) * Gx' * Wꜝx * Hx * uᵧn + convxn + convxn1 + Vx/2 * (fₒn + fₒn1) # Accumulation + Diffusion + Convection (t=0)
    b2 = zeros(n)
    b3 = zeros(n)
    b4 = zeros(n)

    b = [b1; b2; b3; b4]
    return b
end

function solve_NavierStokesUnsteadyMono!(s::VectorSolver, velocity::Velocity, Δt, Tend, bc; method=IterativeSolvers.bicgstabl, kwargs...)
    if s.A === nothing
        error("Solver is not initialized. Call a solver constructor first.")
    end

    n = Int(size(s.A, 1) / 4) # For 1D

    t=0.0

    s.x = method(s.A, s.b; kwargs...)
    push!(s.states, s.x)
    while t<Tend
        t+=Δt
        println("t = ", t)
        # Update the right-hand side
        if typeof(velocity)==Velocity{1}
            s.A = build_navier_stokes_matrix(velocity.operator, velocity.capacities, velocity.ρ, velocity.Re, Δt, s.x, s.states[end-1])
            s.b = build_navier_stokes_rhs(velocity.operator, velocity.capacities, velocity.source, t, velocity.ρ, velocity.Re, Δt, s.x, s.states[end-1])
        elseif typeof(velocity)==Velocity{2}
            println("2D")
        end

        #Apply the boundary conditions
        #BC_border_mono!(s.A, s.b, bc, phase.capacity.mesh)
        
        # Solve the linear system
        s.x = method(s.A, s.b; kwargs...)

        # Update the solution
        push!(s.states, s.x)
    end
end
    


function AdvectionVecUnsteadyMono(velocity, bc, Δt, Tend, x0)
    println("Création du solveur Advection Vectoriel Monophasic Unsteady")
    println("- Monophasic problem")
    println("- Unsteady problem")

    s = VectorSolver(Unsteady, Monophasic, Advection, nothing, nothing, nothing, ConvergenceHistory(), [])

    if typeof(velocity)==Velocity{1}
        println("1D")
    elseif typeof(velocity)==Velocity{2}
        s.A = build_advection_vec_matrix_init(velocity.operator, velocity.capacities, velocity.ρ, Δt, x0)
        s.b = build_advection_vec_rhs_init(velocity.operator, velocity.capacities, velocity.source, Δt, x0)
        push!(s.states, x0)
    end

    # BC_border_mono!(s.A, s.b, bc_b, phase.capacity.mesh)
    return s
end

function build_advection_vec_matrix_init(operator, capacities, ρ, Δt, x0)
    n = prod(operator.size)

    #A = spzeros(Float64, n, n)
    Cx, Cy = operator.C[1], operator.C[2]
    Kx, Ky = operator.K[1], operator.K[2]
    Vx, Vy = operator.V[1], operator.V[2]

    # Compute the blocks
    block1 = Vx/Δt + 1/2 * Cx + 1/4 * Kx
    block2 = 1/4 * Kx
    block3 = I(n)
    block4 = Vy/Δt + 1/2 * Cy + 1/4 * Ky
    block5 = 1/4 * Ky
    block6 = I(n)
    zeroblock = zeros(n, n)

    # Fill the matrix
    A = [block1 block2 zeroblock zeroblock;
        zeroblock block3 zeroblock zeroblock;
        zeroblock zeroblock block4 block5;
        zeroblock zeroblock zeroblock block6]

    return A
end

function build_advection_vec_rhs_init(operator, capacities, source, Δt, x0)
    n = prod(operator.size)
    b = zeros(4n)

    Cx, Cy = operator.C[1], operator.C[2]
    Kx, Ky = operator.K[1], operator.K[2]
    Vx, Vy = operator.V[1], operator.V[2]

    uₒx, uᵧx = x0[1:n], x0[n+1:2n]
    uₒy, uᵧy = x0[2n+1:3n], x0[3n+1:4n]

    b1 = Vx/Δt * uₒx - 1/2 * Cx * uₒx - 1/4 * Kx * uₒx - 1/2 * Kx * uᵧx
    b2 = uᵧx
    b3 = Vy/Δt * uₒy - 1/2 * Cy * uₒy - 1/4 * Ky * uₒy - 1/2 * Ky * uᵧy
    b4 = uᵧy

    b = [b1; b2; b3; b4]
    return b
end

function build_advection_vec_matrix(operator, capacities, ρ, Δt, xn, xn1)
    n = prod(operator.size)
    nx, ny = operator.size

    uₒx, uᵧx = xn[1:n], xn[n+1:2n]
    uₒy, uᵧy = xn[2n+1:3n], xn[3n+1:4n]
    uₒx1, uᵧx1 = xn1[1:n], xn1[n+1:2n]
    uₒy1, uᵧy1 = xn1[2n+1:3n], xn1[3n+1:4n]

    Vx, Vy = operator.V[1], operator.V[2]

    Dx_m, Dy_m = kron(I(ny), ẟ_m(nx)), kron(ẟ_m(ny), I(nx))
    Dx_p, Dy_p = kron(I(ny), δ_p(nx)), kron(δ_p(ny), I(nx))
    Sx_m, Sy_m = kron(I(ny), Σ_m(nx)), kron(Σ_m(ny), I(nx))
    Sx_p, Sy_p = kron(I(ny), Σ_p(nx)), kron(Σ_p(ny), I(nx))

    Cxun = Dx_p * spdiagm(0 => (Sx_m * capacities[1].A[1] * uₒx[1])) * Sx_m 
    Cyun = Dy_p * spdiagm(0 => (Sy_m * capacities[1].A[2] * uₒy[2])) * Sy_m
    Cxvn = Dx_p * spdiagm(0 => (Sx_m * capacities[2].A[1] * uₒ[1])) * Sx_m
    Cyvn = Dy_p * spdiagm(0 => (Sy_m * capacities[2].A[2] * uₒ[2])) * Sy_m
    Cx = Cxu + Cxv
    Cy = Cyu + Cyv

    Kx = spdiagm(0 => Sx_p * Hu' * uᵧ)
    Ky = spdiagm(0 => Sy_p * Hv' * uᵧ)

    # Compute the blocks
    block1 = Vx/Δt 
    block2 = 1/4 * Kx
    block3 = I(n)
    block4 = Vy/Δt + 1/2 * Cy + 3/4 * Ky - 1/4 * Ky
    block5 = 1/4 * Ky
    block6 = I(n)
    zeroblock = zeros(n, n)

    # Fill the matrix
    A = [block1 block2 zeroblock zeroblock;
        zeroblock block3 zeroblock zeroblock;
        zeroblock zeroblock block4 block5;
        zeroblock zeroblock zeroblock block6]

    return A
end