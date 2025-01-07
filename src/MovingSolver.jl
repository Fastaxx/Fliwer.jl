function MovingDiffusionUnsteadyMono(phase::Phase, bc_b::BorderConditions, bc_i::AbstractBoundary, Δt::Float64, Tₑ::Float64, Tᵢ::Vector{Float64})
    println("Création du solveur:")
    println("- Moving problem")
    println("- Monophasic problem")
    println("- Unsteady problem")
    println("- Diffusion problem")
    
    s = Solver(Unsteady, Monophasic, Diffusion, nothing, nothing, nothing, ConvergenceHistory(), [])
    
    s.A = build_mono_unstead_diff_moving_matrix(phase.operator, phase.capacity, phase.Diffusion_coeff, bc_b, bc_i, Δt)
    s.b = build_rhs_mono_unstead_moving_diff(phase.operator, phase.source, phase.capacity, bc_b, bc_i, Tᵢ, Δt, 0.0)

    return s
end

function psip(args::Vararg{T,2}) where {T<:Real}
    if all(iszero, args)
        0.0
    elseif all(!iszero, args)
        0.5
    else
        1.0
    end
end

function psim(args::Vararg{T,2}) where {T<:Real}
    if all(iszero, args)
        0.0
    elseif all(!iszero, args)
        0.5
    else
        0.0
    end
end

function build_mono_unstead_diff_moving_matrix(operator::SpaceTimeOps, capacity::Capacity, D::Float64, bc_b::BorderConditions, bc::AbstractBoundary, Δt::Float64)
    n = prod(operator.size)
    nx, nt = operator.size
    Iₐ, Iᵦ = build_I_bc(operator, bc)
    Iᵧ = build_I_g(operator) # capacite.Γ #
    Id = build_I_D(operator, D)

    Vn_1 = capacity.A[2][1:end÷2, 1:end÷2]
    Vn = capacity.A[2][end÷2+1:end, end÷2+1:end]

    Ψn1 = Diagonal(psip.(Vn,Vn_1))

    W! = operator.Wꜝ[1:nx, 1:nx]
    G = operator.G[1:nx, 1:nx]
    H = operator.H[1:nx, 1:nx]
    Iᵦ = Iᵦ[1:nx, 1:nx]
    Iₐ = Iₐ[1:nx, 1:nx]
    Iᵧ = Iᵧ[1:nx, 1:nx]

    block1 = Vn_1/Δt + G' * W! * G * Ψn1
    block2 = -(Vn_1 - Vn)/Δt + G' * W! * H
    block3 = Iᵦ * H' * W! * G
    block4 = Iᵦ * H' * W! * H + (Iₐ * Iᵧ)

    A = [block1 block2; block3 block4]

    return A
end

function build_rhs_mono_unstead_moving_diff(operator::SpaceTimeOps, f::Function, capacity::Capacity, bc_b::BorderConditions, bc::AbstractBoundary, Tᵢ::Vector{Float64}, Δt::Float64, t::Float64)
    N = prod(operator.size)
    nx,nt = operator.size

    Iᵧ = build_I_g(operator) #capacite.Γ #
    fₒn, fₒn1 = build_source(operator, f, t, capacity), build_source(operator, f, t+Δt, capacity)
    gᵧ = build_g_g(operator, bc, capacity)

    Tₒ, Tᵧ = Tᵢ[1:nx], Tᵢ[nx+1:end]

    Vn_1 = capacity.A[2][1:end÷2, 1:end÷2]
    Vn = capacity.A[2][end÷2+1:end, end÷2+1:end]

    Ψn = Diagonal(psim.(Vn,Vn_1))

    W! = operator.Wꜝ[1:nx, 1:nx]
    G = operator.G[1:nx, 1:nx]
    H = operator.H[1:nx, 1:nx]
    Iᵧ = Iᵧ[1:nx, 1:nx]
    gᵧ = gᵧ[1:nx]

    # Build the right-hand side
    b1 = (Vn/Δt - G' * W! * G * Ψn)*Tₒ - 0.5 * G' * W! * H * Tᵧ # + Vn*(fₒn + fₒn1)/2
    b2 = Iᵧ * gᵧ

    b = [b1; b2]

   return b
end

function solve_MovingDiffusionUnsteadyMono!(s::Solver, phase::Phase, Tᵢ::Vector{Float64}, Δt::Float64, Tₑ::Float64, nt::Int, bc_b::BorderConditions, bc::AbstractBoundary, body::Body, mesh::CartesianMesh, t::Vector{Float64}; method=IterativeSolvers.gmres, kwargs...)
    if s.A === nothing
        error("Solver is not initialized. Call a solver constructor first.")
    end

    println("Solving the problem:")
    println("- Moving problem")
    println("- Monophasic problem")
    println("- Unsteady problem")
    println("- Diffusion problem")

    nx, _ = phase.operator.size
    for i in 2:nt
        println("Time : $(t[i])")
        spaceTimeMesh = CartesianSpaceTimeMesh(mesh, t[i:i+1])
        capacity = Capacity(body, spaceTimeMesh)
        operator = SpaceTimeOps(capacity.A, capacity.B, capacity.V, capacity.W, (nx, 2))
        s.A = build_mono_unstead_diff_moving_matrix(operator, capacity, phase.Diffusion_coeff, bc_b, bc, Δt)
        s.b = build_rhs_mono_unstead_moving_diff(operator, phase.source, capacity, bc_b, bc, Tᵢ, Δt, t[i])
        
        s.x = method(s.A, s.b; kwargs...)

        push!(s.states, s.x)
        @show maximum(s.x)
        Tᵢ = s.x
    end
end