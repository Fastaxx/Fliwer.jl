function MovingDiffusionUnsteadyMono(phase::Phase, bc_b::BorderConditions, bc_i::AbstractBoundary, Δt::Float64, Tₑ::Float64, Tᵢ::Vector{Float64}, scheme::String)
    println("Création du solveur:")
    println("- Moving problem")
    println("- Monophasic problem")
    println("- Unsteady problem")
    println("- Diffusion problem")
    
    s = Solver(Unsteady, Monophasic, Diffusion, nothing, nothing, nothing, ConvergenceHistory(), [])
    
    if scheme == "CN"
        s.A = build_mono_unstead_diff_moving_matrix(phase.operator, phase.capacity, phase.Diffusion_coeff, bc_b, bc_i, Δt, "CN")
        s.b = build_rhs_mono_unstead_moving_diff(phase.operator, phase.source, phase.capacity, bc_b, bc_i, Tᵢ, Δt, 0.0, "CN")
    else 
        s.A = build_mono_unstead_diff_moving_matrix(phase.operator, phase.capacity, phase.Diffusion_coeff, bc_b, bc_i, Δt, "BE")
        s.b = build_rhs_mono_unstead_moving_diff(phase.operator, phase.source, phase.capacity, bc_b, bc_i, Tᵢ, Δt, 0.0, "BE")
    end

    return s
end

function psip_cn(args::Vararg{T,2}) where {T<:Real}
    if all(iszero, args)
        0.0
    elseif all(!iszero, args)
        0.5
    else
        1.0
    end
end

function psim_cn(args::Vararg{T,2}) where {T<:Real}
    if all(iszero, args)
        0.0
    elseif all(!iszero, args)
        0.5
    else
        0.0
    end
end

function psip_be(args::Vararg{T,2}) where {T<:Real}
    if all(iszero, args)
        0.0
    elseif all(!iszero, args)
        1.0
    else
        1.0
    end
end

function psim_be(args::Vararg{T,2}) where {T<:Real}
    if all(iszero, args)
        0.0
    elseif all(!iszero, args)
        0.0
    else
        0.0
    end
end

function build_mono_unstead_diff_moving_matrix(operator::SpaceTimeOps, capacity::Capacity, D::Float64, bc_b::BorderConditions, bc::AbstractBoundary, Δt::Float64, scheme::String)
    n = prod(operator.size)
    nx, nt = operator.size
    Iₐ, Iᵦ = build_I_bc(operator, bc)
    Iᵧ = build_I_g(operator) # capacite.Γ #
    Id = build_I_D(operator, D, capacity)

    Vn_1 = capacity.A[2][1:end÷2, 1:end÷2]
    Vn = capacity.A[2][end÷2+1:end, end÷2+1:end]

    if scheme == "CN"
        psip, psim = psip_cn, psim_cn
    else
        psip, psim = psip_be, psim_be
    end

    Ψn1 = Diagonal(psip.(Vn,Vn_1))

    W! = operator.Wꜝ[1:nx, 1:nx]
    G = operator.G[1:nx, 1:nx]
    H = operator.H[1:nx, 1:nx]
    Iᵦ = Iᵦ[1:nx, 1:nx]
    Iₐ = Iₐ[1:nx, 1:nx]
    Iᵧ = Iᵧ[1:nx, 1:nx]

    block1 = Vn_1 + G' * W! * G * Ψn1
    block2 = -(Vn_1 - Vn) + G' * W! * H * Ψn1 
    block3 = Iᵦ * H' * W! * G 
    block4 = Iᵦ * H' * W! * H + (Iₐ * Iᵧ) 

    A = [block1 block2; block3 block4]
    return A
end

function build_rhs_mono_unstead_moving_diff(operator::SpaceTimeOps, f::Function, capacity::Capacity, bc_b::BorderConditions, bc::AbstractBoundary, Tᵢ::Vector{Float64}, Δt::Float64, t::Float64, scheme::String)
    N = prod(operator.size)
    nx,nt = operator.size

    Iᵧ = build_I_g(operator) #capacite.Γ #
    fₒn, fₒn1 = build_source(operator, f, t, capacity), build_source(operator, f, t+Δt, capacity)
    gᵧ = build_g_g(operator, bc, capacity)

    Tₒ, Tᵧ = Tᵢ[1:nx], Tᵢ[nx+1:end]

    Vn_1 = capacity.A[2][1:end÷2, 1:end÷2]
    Vn = capacity.A[2][end÷2+1:end, end÷2+1:end]

    if scheme == "CN"
        psip, psim = psip_cn, psim_cn
    else
        psip, psim = psip_be, psim_be
    end

    Ψn = Diagonal(psim.(Vn,Vn_1))

    W! = operator.Wꜝ[1:nx, 1:nx]
    G = operator.G[1:nx, 1:nx]
    H = operator.H[1:nx, 1:nx]
    V = operator.V[1:nx, 1:nx]
    Iᵧ = Iᵧ[1:nx, 1:nx]
    gᵧ = gᵧ[1:nx]
    fₒn, fₒn1 = fₒn[1:nx], fₒn1[1:nx]

    # Build the right-hand side
    if scheme == "CN"
        b1 = (Vn - G' * W! * G * Ψn)*Tₒ - 0.5 * G' * W! * H * Tᵧ + 0.5 * V * (fₒn + fₒn1)
    else
        b1 = (Vn)*Tₒ + V * (fₒn1)
    end
    b2 = Iᵧ * gᵧ

    b = [b1; b2]

   return b
end

function solve_MovingDiffusionUnsteadyMono!(
    s::Solver,
    phase::Phase,
    Tᵢ::Vector{Float64},
    Δt::Float64,
    Tₑ::Float64,
    nt::Int,
    bc_b::BorderConditions,
    bc::AbstractBoundary,
    body::Body,
    mesh::CartesianMesh,
    t::Vector{Float64},
    scheme::String;
    method = IterativeSolvers.gmres,
    kwargs...
)

    if s.A === nothing
        error("Solver is not initialized. Call a solver constructor first.")
    end

    println("Solving the problem:")
    println("- Moving problem")
    println("- Monophasic problem")
    println("- Unsteady problem")
    println("- Diffusion problem")

    nx, _ = phase.operator.size
    cond_log = Float64[]
    minV_log = Float64[]
    maxV_log = Float64[]
    minW_log = Float64[]
    maxW_log = Float64[]

    for i in 2:nt
        println("Time : $(t[i])")
        spaceTimeMesh = CartesianSpaceTimeMesh(mesh, t[i:i+1];tag=mesh.tag)
        capacity = Capacity(body, spaceTimeMesh)
        operator = SpaceTimeOps(
            capacity.A, capacity.B,
            capacity.V, capacity.W,
            (nx, 2)
        )

        s.A = build_mono_unstead_diff_moving_matrix(
            operator,
            capacity,
            phase.Diffusion_coeff,
            bc_b,
            bc,
            Δt,
            scheme
        )
        s.b = build_rhs_mono_unstead_moving_diff(
            operator,
            phase.source,
            capacity,
            bc_b,
            bc,
            Tᵢ,
            Δt,
            t[i],
            scheme
        )
        BC_border_mono!(s.A, s.b, bc_b, capacity.mesh)



        # Solve system
        if method == \
            A_reduced, b_reduced, rows_idx, cols_idx = remove_zero_rows_cols!(s.A, s.b)
            # Compute condition number
            cnum = cond(Array(A_reduced), 2)
            push!(cond_log, cnum)
            push!(minV_log, minimum(x for x in capacity.V if x != 0))
            push!(maxV_log, maximum(capacity.V))
            push!(minW_log, minimum(x for x in capacity.W[1] if x != 0))
            push!(maxW_log, maximum(capacity.W[1]))
            x_reduced = A_reduced \ b_reduced
            s.x = zeros(size(s.A, 1))
            s.x[cols_idx] .= x_reduced
        else
            log = get(kwargs, :log, false)
            if log
                s.x, s.ch = method(s.A, s.b; kwargs...)
            else
                s.x = method(s.A, s.b; kwargs...)
            end
        end

        push!(s.states, s.x)
        @show maximum(abs.(s.x))
        Tᵢ = s.x
    end

    # Store condition numbers & min/max in a file
    open("condition_numbers.txt", "w") do f
        for j in 1:length(cond_log)
            println(f, join((
                cond_log[j],
                minV_log[j],
                maxV_log[j],
                minW_log[j],
                maxW_log[j]
            ), " "))
        end
    end
end