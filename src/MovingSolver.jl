# Moving - Diffusion - Unsteady - Monophasic
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

"""
function psip_cn(args::Vararg{T,2}) where {T<:Real}
    if all(iszero, args)
        0.0
    elseif all(!iszero, args)
        0.5
    else
        0.5
    end
end
"""

function psip_cn(args::Vararg{T,2}) where {T<:Real}
    if all(iszero, args)
        0.0
    elseif all(!iszero, args)
        0.5
    elseif iszero(args[1]) && !iszero(args[2])
        0.5
    elseif !iszero(args[1]) && iszero(args[2])
        1.0
    else
        0.0
    end        
end

function psim_cn(args::Vararg{T,2}) where {T<:Real}
    if all(iszero, args)
        0.0
    elseif all(!iszero, args)
        0.5
    elseif iszero(args[1]) && !iszero(args[2]) # Fresh
        0.5
    elseif !iszero(args[1]) && iszero(args[2]) # Dead
        0.0
    else
        0.0
    end
end


"""
function psim_cn(args::Vararg{T,2}) where {T<:Real}
    if all(iszero, args)
        0.0
    elseif all(!iszero, args)
        0.5
    else
        0.5
    end
end
"""

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

        # CFL log
        Vn_1 = capacity.A[2][1:end÷2, 1:end÷2]
        Vn = capacity.A[2][end÷2+1:end, end÷2+1:end]
        cfln1 = 0.2 * Δt./(maximum(Vn_1))
        cfln = 0.2 * Δt./(maximum(Vn))

        println("CFL number at time n+1 : ", cfln1)
        println("CFL number at time n : ", cfln)    

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


function MovingDiffusionUnsteadyMono2(phase::Phase, bc_b::BorderConditions, bc_i::AbstractBoundary, Δt::Float64, Tₑ::Float64, Tᵢ::Vector{Float64}, scheme::String)
    println("Création du solveur:")
    println("- Moving problem")
    println("- Monophasic problem")
    println("- Unsteady problem")
    println("- Diffusion problem")
    
    s = Solver(Unsteady, Monophasic, Diffusion, nothing, nothing, nothing, ConvergenceHistory(), [])
    
    if scheme == "CN"
        s.A = build_mono_unstead_diff_moving_matrix2(phase.operator, phase.capacity, phase.Diffusion_coeff, bc_b, bc_i, Δt, "CN")
        s.b = build_rhs_mono_unstead_moving_diff2(phase.operator, phase.source, phase.capacity, bc_b, bc_i, Tᵢ, Δt, 0.0, "CN")
    else 
        s.A = build_mono_unstead_diff_moving_matrix2(phase.operator, phase.capacity, phase.Diffusion_coeff, bc_b, bc_i, Δt, "BE")
        s.b = build_rhs_mono_unstead_moving_diff2(phase.operator, phase.source, phase.capacity, bc_b, bc_i, Tᵢ, Δt, 0.0, "BE")
    end

    return s
end

function build_mono_unstead_diff_moving_matrix2(operator::SpaceTimeOps, capacity::Capacity, D::Float64, bc_b::BorderConditions, bc::AbstractBoundary, Δt::Float64, scheme::String)
    nx, ny, nt = operator.size
    n = nx*ny
    Iₐ, Iᵦ = build_I_bc(operator, bc)
    Iᵧ = build_I_g(operator) # capacite.Γ #
    Id = build_I_D(operator, D, capacity)

    Vn_1 = capacity.A[3][1:end÷2, 1:end÷2]
    Vn = capacity.A[3][end÷2+1:end, end÷2+1:end]
    if scheme == "CN"
        psip, psim = psip_cn, psim_cn
    else
        psip, psim = psip_be, psim_be
    end

    Ψn1 = Diagonal(psip.(Vn,Vn_1))

    W! = operator.Wꜝ[1:n, 1:n]
    G = operator.G[1:n, 1:n]
    H = operator.H[1:n, 1:n]
    Iᵦ = Iᵦ[1:n, 1:n]
    Iₐ = Iₐ[1:n, 1:n]
    Iᵧ = Iᵧ[1:n, 1:n]

    block1 = Vn_1 + G' * W! * G * Ψn1
    block2 = -(Vn_1 - Vn) + G' * W! * H * Ψn1
    block3 = Iᵦ * H' * W! * G 
    block4 = Iᵦ * H' * W! * H + (Iₐ * Iᵧ) 

    A = [block1 block2; block3 block4]
    return A
end

function build_rhs_mono_unstead_moving_diff2(operator::SpaceTimeOps, f::Function, capacity::Capacity, bc_b::BorderConditions, bc::AbstractBoundary, Tᵢ::Vector{Float64}, Δt::Float64, t::Float64, scheme::String)
    nx, ny, nt = operator.size
    n = nx*ny
    Iᵧ = build_I_g(operator) #capacite.Γ #
    fₒn, fₒn1 = build_source(operator, f, t, capacity), build_source(operator, f, t+Δt, capacity)
    gᵧ = build_g_g(operator, bc, capacity)

    Tₒ, Tᵧ = Tᵢ[1:n], Tᵢ[n+1:end]
    Vn_1 = capacity.A[3][1:end÷2, 1:end÷2]
    Vn = capacity.A[3][end÷2+1:end, end÷2+1:end]

    if scheme == "CN"
        psip, psim = psip_cn, psim_cn
    else
        psip, psim = psip_be, psim_be
    end

    Ψn = Diagonal(psim.(Vn,Vn_1))

    W! = operator.Wꜝ[1:n, 1:n]
    G = operator.G[1:n, 1:n]
    H = operator.H[1:n, 1:n]
    V = operator.V[1:n, 1:n]
    Iᵧ = Iᵧ[1:n, 1:n]
    gᵧ = gᵧ[1:n]
    fₒn, fₒn1 = fₒn[1:n], fₒn1[1:n]

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

function solve_MovingDiffusionUnsteadyMono2!(
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

    nx, ny, _ = phase.operator.size
    cond_log = Float64[]
    minV_log = Float64[]
    maxV_log = Float64[]
    minW_log = Float64[]
    maxW_log = Float64[]
    maxT_log = Float64[]

    # Solve for the initial condition
    BC_border_mono!(s.A, s.b, bc_b, phase.capacity.mesh)

    # Solve system
    if method == \
        A_reduced, b_reduced, rows_idx, cols_idx = remove_zero_rows_cols!(s.A, s.b)
        # Compute condition number
        #cnum = cond(Array(A_reduced), 2)
        cnum = 0.0
        push!(cond_log, cnum)
        push!(minV_log, minimum(x for x in phase.capacity.V if x != 0))
        push!(maxV_log, maximum(phase.capacity.V))
        push!(minW_log, minimum(x for x in phase.capacity.W[1] if x != 0))
        push!(maxW_log, maximum(phase.capacity.W[1]))
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
    push!(maxT_log, maximum(abs.(s.x)))
    Tᵢ = s.x

    # Solve for the second time step
    println("Time : $(t[2])")
    spaceTimeMesh = CartesianSpaceTimeMesh(mesh, t[2:3];tag=mesh.tag)
    capacity = Capacity(body, spaceTimeMesh)
    operator = SpaceTimeOps(
        capacity.A, capacity.B,
        capacity.V, capacity.W,
        (nx, ny, 2)
    )

    s.A = build_mono_unstead_diff_moving_matrix2(
        operator,
        capacity,
        phase.Diffusion_coeff,
        bc_b,
        bc,
        Δt,
        "BE"
    )

    s.b = build_rhs_mono_unstead_moving_diff2(
        operator,
        phase.source,
        capacity,
        bc_b,
        bc,
        Tᵢ,
        Δt,
        t[2],
        "BE"
    )

    BC_border_mono!(s.A, s.b, bc_b, capacity.mesh)

    # Solve system
    if method == \
        A_reduced, b_reduced, rows_idx, cols_idx = remove_zero_rows_cols!(s.A, s.b)
        # Compute condition number
        #cnum = cond(Array(A_reduced), 2)
        cnum = 0.0
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
    push!(maxT_log, maximum(abs.(s.x)))
    Tᵢ = s.x

    for i in 3:nt
        println("Time : $(t[i])")
        spaceTimeMesh = CartesianSpaceTimeMesh(mesh, t[i:i+1];tag=mesh.tag)
        capacity = Capacity(body, spaceTimeMesh)
        operator = SpaceTimeOps(
            capacity.A, capacity.B,
            capacity.V, capacity.W,
            (nx, ny, 2)
        )

        s.A = build_mono_unstead_diff_moving_matrix2(
            operator,
            capacity,
            phase.Diffusion_coeff,
            bc_b,
            bc,
            Δt,
            scheme
        )
        s.b = build_rhs_mono_unstead_moving_diff2(
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

        # CFL log
        Vn_1 = capacity.A[2][1:end÷2, 1:end÷2]
        Vn = capacity.A[2][end÷2+1:end, end÷2+1:end]
        cflmax = 0.02 * Δt./(maximum(Vn_1))
        cflmin = 0.02 * Δt./(minimum(Vn_1))

        println("CFL number max : ", cflmax)
        println("CFL number min : ", cflmin)   

        # Solve system
        if method == \
            A_reduced, b_reduced, rows_idx, cols_idx = remove_zero_rows_cols!(s.A, s.b)
            # Compute condition number
            #cnum = cond(Array(A_reduced), 2)
            cnum = 0.0
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
        push!(maxT_log, maximum(abs.(s.x)))
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

    # Store maxT in a file
    open("max_T_log.txt", "w") do f
        for j in 1:length(maxT_log)
            println(f, maxT_log[j])
        end
    end
end










# Moving - Diffusion - Unsteady - Diphasic
function MovingDiffusionUnsteadyDiph(phase1::Phase, phase2::Phase, bc_b::BorderConditions, ic::InterfaceConditions, Δt::Float64, Tₑ::Float64, Tᵢ::Vector{Float64}, scheme::String)
    println("Création du solveur:")
    println("- Moving problem")
    println("- Diphasic problem")
    println("- Unsteady problem")
    println("- Diffusion problem")
    
    s = Solver(Unsteady, Diphasic, Diffusion, nothing, nothing, nothing, ConvergenceHistory(), [])
    
    if scheme == "CN"
        s.A = build_diph_unstead_diff_moving_matrix(phase1.operator, phase2.operator, phase1.capacity, phase2.capacity, phase1.Diffusion_coeff, phase2.Diffusion_coeff, bc_b, ic, Δt, "CN")
        s.b = build_rhs_diph_unstead_moving_diff(phase1.operator, phase2.operator, phase1.source, phase2.source, phase1.capacity, phase2.capacity, bc_b, ic, Tᵢ, Δt, 0.0, "CN")
    else 
        s.A = build_diph_unstead_diff_moving_matrix(phase1.operator, phase2.operator, phase1.capacity, phase2.capacity, phase1.Diffusion_coeff, phase2.Diffusion_coeff, bc_b, ic, Δt, "BE")
        s.b = build_rhs_diph_unstead_moving_diff(phase1.operator, phase2.operator, phase1.source, phase2.source, phase1.capacity, phase2.capacity, bc_b, ic, Tᵢ, Δt, 0.0, "BE")
    end

    return s
end

function build_diph_unstead_diff_moving_matrix(operator1::SpaceTimeOps, operator2::SpaceTimeOps, capacite1::Capacity, capacite2::Capacity, D1::Float64, D2::Float64, bc_b::BorderConditions, ic::InterfaceConditions, Δt::Float64, scheme::String)
    n1, n2 = prod(operator1.size), prod(operator2.size)
    nx1, nt1 = operator1.size
    nx2, nt2 = operator2.size

    jump, flux = ic.scalar, ic.flux
    Iₐ1, Iₐ2 = jump.α₁ * I(n1), jump.α₂ * I(n2)
    Iᵦ1, Iᵦ2 = flux.β₁ * I(n2), flux.β₂ * I(n2)
    Id1, Id2 = build_I_D(operator1, D1, capacite1), build_I_D(operator2, D2, capacite2)

    Vn1_1 = capacite1.A[2][1:end÷2, 1:end÷2]
    Vn1 = capacite1.A[2][end÷2+1:end, end÷2+1:end]
    Vn2_1 = capacite2.A[2][1:end÷2, 1:end÷2]
    Vn2 = capacite2.A[2][end÷2+1:end, end÷2+1:end]

    if scheme == "CN"
        psip, psim = psip_cn, psim_cn
    else
        psip, psim = psip_be, psim_be
    end

    Ψn1 = Diagonal(psip.(Vn1,Vn1_1))
    Ψn2 = Diagonal(psip.(Vn2,Vn2_1))

    W!1 = operator1.Wꜝ[1:nx1, 1:nx1]
    G1 = operator1.G[1:nx1, 1:nx1]
    H1 = operator1.H[1:nx1, 1:nx1]
    W!2 = operator2.Wꜝ[1:nx2, 1:nx2]
    G2 = operator2.G[1:nx2, 1:nx2]
    H2 = operator2.H[1:nx2, 1:nx2]
    Iᵦ1 = Iᵦ1[1:nx1, 1:nx1]
    Iᵦ2 = Iᵦ2[1:nx2, 1:nx2]
    Iₐ1 = Iₐ1[1:nx1, 1:nx1]
    Iₐ2 = Iₐ2[1:nx2, 1:nx2]

    block1 = Vn1_1 + G1' * W!1 * G1 * Ψn1
    block2 = -(Vn1_1 - Vn1) + G1' * W!1 * H1 * Ψn1
    block3 = Vn2_1 + G2' * W!2 * G2 * Ψn2
    block4 = -(Vn2_1 - Vn2) + G2' * W!2 * H2 * Ψn2
    block5 = Iᵦ1 * H1' * W!1 * G1
    block6 = Iᵦ1 * H1' * W!1 * H1
    block7 = Iᵦ2 * H2' * W!2 * G2
    block8 = Iᵦ2 * H2' * W!2 * H2

    n=nx1
    # Preallocate the sparse matrix
    A = spzeros(Float64, 4n, 4n)
    
    # Assign blocks to the matrix
    A[1:n, 1:n] = block1
    A[1:n, n+1:2n] = block2
    A[1:n, 2n+1:3n] = spzeros(n, n)
    A[1:n, 3n+1:4n] = spzeros(n, n)

    A[n+1:2n, 1:n] = spzeros(n, n)
    A[n+1:2n, n+1:2n] = Iₐ1
    A[n+1:2n, 2n+1:3n] = spzeros(n, n)
    A[n+1:2n, 3n+1:4n] = -Iₐ2

    A[2n+1:3n, 1:n] = spzeros(n, n)
    A[2n+1:3n, n+1:2n] = spzeros(n, n)
    A[2n+1:3n, 2n+1:3n] = block3
    A[2n+1:3n, 3n+1:4n] = block4

    A[3n+1:4n, 1:n] = block5
    A[3n+1:4n, n+1:2n] = block6
    A[3n+1:4n, 2n+1:3n] = block7
    A[3n+1:4n, 3n+1:4n] = block8


    return A
end

function build_rhs_diph_unstead_moving_diff(operator1::SpaceTimeOps, operator2::SpaceTimeOps, f1::Function, f2::Function, capacite1::Capacity, capacite2::Capacity, bc_b::BorderConditions, ic::InterfaceConditions, Tᵢ::Vector{Float64}, Δt::Float64, t::Float64, scheme::String)
    n1, n2 = prod(operator1.size), prod(operator2.size)
    nx1, nt1 = operator1.size
    nx2, nt2 = operator2.size

    jump, flux = ic.scalar, ic.flux
    Iₐ1, Iₐ2 = jump.α₁ * I(n1), jump.α₂ * I(n2)
    Iᵦ1, Iᵦ2 = flux.β₁ * I(n2), flux.β₂ * I(n2)

    f1ₒn, f1ₒn1 = build_source(operator1, f1, t, capacite1), build_source(operator1, f1, t+Δt, capacite1)
    f2ₒn, f2ₒn1 = build_source(operator2, f2, t, capacite2), build_source(operator2, f2, t+Δt, capacite2)

    Iᵧ1, Iᵧ2 = build_I_g(operator1), build_I_g(operator2)
    gᵧ, hᵧ = build_g_g(operator1, jump, capacite1), build_g_g(operator2, flux, capacite2)

    Vn1_1 = capacite1.A[2][1:end÷2, 1:end÷2]
    Vn1 = capacite1.A[2][end÷2+1:end, end÷2+1:end]
    Vn2_1 = capacite2.A[2][1:end÷2, 1:end÷2]
    Vn2 = capacite2.A[2][end÷2+1:end, end÷2+1:end]

    if scheme == "CN"
        psip, psim = psip_cn, psim_cn
    else
        psip, psim = psip_be, psim_be
    end

    Ψn1 = Diagonal(psim.(Vn1,Vn1_1))
    Ψn2 = Diagonal(psim.(Vn2,Vn2_1))

    W!1 = operator1.Wꜝ[1:nx1, 1:nx1]
    G1 = operator1.G[1:nx1, 1:nx1]
    H1 = operator1.H[1:nx1, 1:nx1]
    W!2 = operator2.Wꜝ[1:nx2, 1:nx2]
    G2 = operator2.G[1:nx2, 1:nx2]
    H2 = operator2.H[1:nx2, 1:nx2]
    V1 = operator1.V[1:nx1, 1:nx1]
    V2 = operator2.V[1:nx2, 1:nx2]
    Iᵦ1 = Iᵦ1[1:nx1, 1:nx1]
    Iᵦ2 = Iᵦ2[1:nx2, 1:nx2]
    Iₐ1 = Iₐ1[1:nx1, 1:nx1]
    Iₐ2 = Iₐ2[1:nx2, 1:nx2]

    Tₒ1, Tᵧ1 = Tᵢ[1:nx1], Tᵢ[nx1+1:2nx1]
    Tₒ2, Tᵧ2 = Tᵢ[2nx1+1:2nx1+nx2], Tᵢ[2nx1+nx2+1:end]

    f1ₒn, f1ₒn1 = f1ₒn[1:nx1], f1ₒn1[1:nx1]
    f2ₒn, f2ₒn1 = f2ₒn[1:nx2], f2ₒn1[1:nx2]
    gᵧ, hᵧ = gᵧ[1:nx1], hᵧ[1:nx2]
    Iᵧ1, Iᵧ2 = Iᵧ1[1:nx1, 1:nx1], Iᵧ2[1:nx2, 1:nx2]

    # Build the right-hand side
    if scheme == "CN"
        b1 = (Vn1 - G1' * W!1 * G1 * Ψn1)*Tₒ1 - 0.5 * G1' * W!1 * H1 * Tᵧ1 + 0.5 * V1 * (f1ₒn + f1ₒn1)
        b3 = (Vn2 - G2' * W!2 * G2 * Ψn2)*Tₒ2 - 0.5 * G2' * W!2 * H2 * Tᵧ2 + 0.5 * V2 * (f2ₒn + f2ₒn1)
    else
        b1 = (Vn1)*Tₒ1 + V1 * (f1ₒn1)
        b3 = (Vn2)*Tₒ2 + V2 * (f2ₒn1)
    end
    b2 = gᵧ
    b4 = Iᵧ2 * hᵧ

    b = [b1; b2; b3; b4]

    return b
end

function solve_MovingDiffusionUnsteadyDiph!(
    s::Solver,
    phase1::Phase,
    phase2::Phase,
    Tᵢ::Vector{Float64},
    Δt::Float64,
    Tₑ::Float64,
    nt::Int,
    bc_b::BorderConditions,
    ic::InterfaceConditions,
    body::Body,
    body_c::Body,
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
    println("- Diphasic problem")
    println("- Unsteady problem")
    println("- Diffusion problem")

    nx1, _ = phase1.operator.size
    nx2, _ = phase2.operator.size
    cond_log = Float64[]
    minV_log = Float64[]
    maxV_log = Float64[]
    minW_log = Float64[]
    maxW_log = Float64[]
    maxT_log = Float64[]

    # Solve for the initial condition
    BC_border_diph!(s.A, s.b, bc_b, phase1.capacity.mesh)

    # Solve system
    if method == \
        A_reduced, b_reduced, rows_idx, cols_idx = remove_zero_rows_cols!(s.A, s.b)
        # Compute condition number
        #cnum = cond(Array(A_reduced), 2)
        cnum = 0.0
        push!(cond_log, cnum)
        push!(minV_log, minimum(x for x in phase1.capacity.V if x != 0))
        push!(maxV_log, maximum(phase1.capacity.V))
        push!(minW_log, minimum(x for x in phase1.capacity.W[1] if x != 0))
        push!(maxW_log, maximum(phase1.capacity.W[1]))
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
    push!(maxT_log, maximum(abs.(s.x)))
    Tᵢ = s.x

    for i in 2:nt
        println("Time : $(t[i])")
        spaceTimeMesh = CartesianSpaceTimeMesh(mesh, t[i:i+1];tag=mesh.tag)
        capacity1 = Capacity(body, spaceTimeMesh)
        operator1 = SpaceTimeOps(
            capacity1.A, capacity1.B,
            capacity1.V, capacity1.W,
            (nx1, 2)
        )

        capacity2 = Capacity(body_c, spaceTimeMesh)
        operator2 = SpaceTimeOps(
            capacity2.A, capacity2.B,
            capacity2.V, capacity2.W,
            (nx2, 2)
        )

        s.A = build_diph_unstead_diff_moving_matrix(
            operator1,
            operator2,
            capacity1,
            capacity2,
            phase1.Diffusion_coeff,
            phase2.Diffusion_coeff,
            bc_b,
            ic,
            Δt,
            scheme
        )
        s.b = build_rhs_diph_unstead_moving_diff(
            operator1,
            operator2,
            phase1.source,
            phase2.source,
            capacity1,
            capacity2,
            bc_b,
            ic,
            Tᵢ,
            Δt,
            t[i],
            scheme
        )
        BC_border_diph!(s.A, s.b, bc_b, mesh)

        # Solve system
        if method == \
            A_reduced, b_reduced, rows_idx, cols_idx = remove_zero_rows_cols!(s.A, s.b)
            # Compute condition number
            #cnum = cond(Array(A_reduced), 2)
            cnum = 0.0
            push!(cond_log, cnum)
            push!(minV_log, minimum(x for x in capacity1.V if x != 0))
            push!(maxV_log, maximum(capacity1.V))
            push!(minW_log, minimum(x for x in capacity1.W[1] if x != 0))
            push!(maxW_log, maximum(capacity1.W[1]))
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
        push!(maxT_log, maximum(abs.(s.x)))
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

    # Store maxT in a file
    open("max_T_log.txt", "w") do f
        for j in 1:length(maxT_log)
            println(f, maxT_log[j])
        end
    end
end