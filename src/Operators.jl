"""
    abstract type AbstractOperators

An abstract type representing a collection of operators.
"""
abstract type AbstractOperators end

"""
    ∇(operator::AbstractOperators, p::Vector{Float64})

Compute the gradient of a scalar field.
"""
function ∇(operator::AbstractOperators, p::Vector{Float64})
    ∇ = operator.Wꜝ * (operator.G * p[1:div(end,2)] + operator.H * p[div(end,2)+1:end])
    return ∇
end

"""
    ∇_(operator::AbstractOperators, qω::Vector{Float64}, qγ::Vector{Float64})

Compute the divergence of a vector field.
"""
function ∇_(operator::AbstractOperators, qω::Vector{Float64}, qγ::Vector{Float64})
    GT = operator.G'
    HT = operator.H'
    return -(GT + HT)*qω + HT * qγ
end

"""
    struct DiffusionOps{N} <: AbstractOperators where N

Struct representing diffusion operators.

# Fields
- `G::SparseMatrixCSC{Float64, Int}`: Matrix representing the diffusion operator G.
- `H::SparseMatrixCSC{Float64, Int}`: Matrix representing the diffusion operator H.
- `Wꜝ::SparseMatrixCSC{Float64, Int}`: Matrix representing the diffusion operator Wꜝ.
- `V::SparseMatrixCSC{Float64, Int}`: Matrix representing the diffusion operator V.
- `size::NTuple{N, Int}`: Tuple representing the size of the diffusion operators.

"""
struct DiffusionOps{N} <: AbstractOperators where N
    G::SparseMatrixCSC{Float64, Int}
    H::SparseMatrixCSC{Float64, Int}
    Wꜝ::SparseMatrixCSC{Float64, Int}
    V::SparseMatrixCSC{Float64, Int}
    size::NTuple{N, Int}
end

# Elementary operators
function ẟ_m(n::Int, periodicity::Bool=false) D = spdiagm(0 => ones(n), -1 => -ones(n-1)); D[n, n] = 0.0; if periodicity; D[1, n-1] = -1.0; D[n, 1] = 1.0; end; D end
function δ_p(n::Int, periodicity::Bool=false) D = spdiagm(0 => -ones(n), 1 => ones(n-1)); D[n, n] = 0.0; if periodicity; D[1, n-1] = -1.0; D[n, 1] = 1.0; end; D end
function Σ_m(n::Int, periodicity::Bool=false) D = 0.5 * spdiagm(0 => ones(n), -1 => ones(n-1)); D[n, n] = 0.0; if periodicity; D[1, n-1] = 0.5; D[n, 1] = 0.5; end; D end
function Σ_p(n::Int, periodicity::Bool=false) D = 0.5 * spdiagm(0 => ones(n), 1 => ones(n-1)); D[n, n] = 0.0; if periodicity; D[1, n-1] = 0.5; D[n, 1] = 0.5; end; D end
function I(n::Int) spdiagm(0 => ones(n)) end

"""
    DiffusionOps(A, B, V, W, size)

Constructs the diffusion operators for a given system.

# Arguments
- `A`: Array of A matrices for each dimension.
- `B`: Array of B matrices for each dimension.
- `V`: V matrix for the system.
- `W`: W matrix for the system.
- `size`: Array of sizes for each dimension.

# Returns
- `DiffusionOps`: Diffusion operators for the system.
"""
function DiffusionOps(A, B, V, W, size)
    N = length(size)
    if N == 1
        nx = size[1]
        G = ẟ_m(nx) * B[1] # Gérer le periodicity
        H = A[1]*ẟ_m(nx) - ẟ_m(nx)*B[1]
        diagW = diag(blockdiag(W[1]))
        new_diagW = [val != 0 ? 1.0 / val : 1.0 for val in diagW]
        Wꜝ = spdiagm(0 => new_diagW)
    elseif N == 2
        nx, ny = size[1], size[2]
        Dx_m = kron(I(ny), ẟ_m(nx))
        Dy_m = kron(ẟ_m(ny), I(nx))
        G = [Dx_m * B[1]; Dy_m * B[2]]
        H = [A[1]*Dx_m - Dx_m*B[1]; A[2]*Dy_m - Dy_m*B[2]]
        diagW = diag(blockdiag(W[1], W[2]))
        new_diagW = [val != 0 ? 1.0 / val : 1.0 for val in diagW]
        Wꜝ = spdiagm(0 => new_diagW)
    elseif N == 3
        nx, ny, nz = size[1], size[2], size[3]
        Dx_m = kron(I(nz), kron(I(ny), ẟ_m(nx)))
        Dy_m = kron(I(nz), kron(ẟ_m(ny), I(nx)))
        Dz_m = kron(ẟ_m(nz), kron(I(ny), I(nx)))
        G = [Dx_m * B[1]; Dy_m * B[2]; Dz_m * B[3]]
        H = [A[1]*Dx_m - Dx_m*B[1]; A[2]*Dy_m - Dy_m*B[2]; A[3]*Dz_m - Dz_m*B[3]]
        diagW = diag(blockdiag(W[1], W[2], W[3]))
        new_diagW = [val != 0 ? 1.0 / val : 1.0 for val in diagW]
        Wꜝ = spdiagm(0 => new_diagW)
    end

    return DiffusionOps{N}(G, H, Wꜝ, V, size)
end

"""
    struct ConvectionOps{N} <: AbstractOperators where N

Struct representing a collection of convection operators.

# Fields
- `C`: A tuple of N sparse matrices representing the C operators.
- `K`: A tuple of N sparse matrices representing the K operators.
- `size`: A tuple of N integers representing the size of each operator.

"""
struct ConvectionOps{N} <: AbstractOperators where N
    C :: NTuple{N, SparseMatrixCSC{Float64, Int}}
    K :: NTuple{N, SparseMatrixCSC{Float64, Int}}
    G::SparseMatrixCSC{Float64, Int}
    H::SparseMatrixCSC{Float64, Int}
    Wꜝ::SparseMatrixCSC{Float64, Int}
    V::SparseMatrixCSC{Float64, Int}
    size :: NTuple{N, Int}
end

"""
    ConvectionOps(C, K, size)

Constructs the convection operators for a given system.
# Arguments
- `A`: Array of A matrices for each dimension.
- `B`: Array of B matrices for each dimension.
- `V`: V matrix for the system.
- `W`: W matrix for the system.
- `size`: Array of sizes for each dimension.
- `uₒ`: Array of uₒ values for each dimension.
- `uᵧ`: Array of uᵧ values for each dimension.

# Returns
- `ConvectionOps`: Convection operators for the system.
"""
function ConvectionOps(A, B, V, W, size, uₒ, uᵧ)
    N = length(size)
    if N == 1
        nx = size[1]
        Cx = δ_p(nx) * spdiagm(0 => (Σ_m(nx) * A[1] * uₒ[1])) * Σ_m(nx)
        G = ẟ_m(nx) * B[1] # Gérer le periodicity
        H = A[1]*ẟ_m(nx) - ẟ_m(nx)*B[1]
        Kx = spdiagm(0 => Σ_p(nx) * H' * uᵧ)
        diagW = diag(blockdiag(W[1]))
        new_diagW = [val != 0 ? 1.0 / val : 1.0 for val in diagW]
        Wꜝ = spdiagm(0 => new_diagW)
        return ConvectionOps{N}((Cx,), (Kx,), G, H, Wꜝ, V, size)
    elseif N == 2
        nx, ny = size[1], size[2]
        Dx_m, Dy_m = kron(I(ny), ẟ_m(nx)), kron(ẟ_m(ny), I(nx))
        Dx_p, Dy_p = kron(I(ny), δ_p(nx)), kron(δ_p(ny), I(nx))
        Sx_p, Sy_p = kron(I(ny), Σ_p(nx)), kron(Σ_p(ny), I(nx))
        Sx_m, Sy_m = kron(I(ny), Σ_m(nx)), kron(Σ_m(ny), I(nx))
        G = [Dx_m * B[1]; Dy_m * B[2]]
        Cx = Dx_p * spdiagm(0 => (Sx_m * A[1] * uₒ[1])) * Sx_m
        Cy = Dy_p * spdiagm(0 => (Sy_m * A[2] * uₒ[2])) * Sy_m
        H = [A[1]*Dx_m - Dx_m*B[1]; A[2]*Dy_m - Dy_m*B[2]]
        Kx = spdiagm(0 => Sx_p * H' * uᵧ)
        Ky = spdiagm(0 => Sy_p * H' * uᵧ)
        diagW = diag(blockdiag(W[1], W[2]))
        new_diagW = [val != 0 ? 1.0 / val : 1.0 for val in diagW]
        Wꜝ = spdiagm(0 => new_diagW)
        return ConvectionOps{N}((Cx, Cy), (Kx, Ky), G, H, Wꜝ, V, size)
    elseif N == 3
        nx, ny, nz = size[1], size[2], size[3]
        Dx_m, Dy_m, Dz_m = kron(I(nz), kron(I(ny), ẟ_m(nx))), kron(I(nz), kron(ẟ_m(ny), I(nx))), kron(ẟ_m(nz), kron(I(ny), I(nx)))
        Sx_m, Sy_m, Sz_m = kron(I(nz), kron(I(ny), Σ_m(nx))), kron(I(nz), kron(Σ_m(ny), I(nx))), kron(Σ_m(nz), kron(I(ny), I(nx)))
        Sx_p, Sy_p, Sz_p = kron(I(nz), kron(I(ny), Σ_p(nx))), kron(I(nz), kron(Σ_p(ny), I(nx))), kron(Σ_p(nz), kron(I(ny), I(nx)))
        G = [Dx_m * B[1]; Dy_m * B[2]; Dz_m * B[3]]
        Cx = Dx_p * spdiagm(0 => (Sx_m * A[1] * uₒ[1])) * Sx_m
        Cy = Dy_p * spdiagm(0 => (Sy_m * A[2] * uₒ[2])) * Sy_m
        Cz = Dz_p * spdiagm(0 => (Sz_m * A[3] * uₒ[3])) * Sz_m
        H = [A[1]*Dx_m - Dx_m*B[1]; A[2]*Dy_m - Dy_m*B[2]; A[3]*Dz_m - Dz_m*B[3]]
        Kx = spdiagm(0 => Sx_p * H' * uᵧ)
        Ky = spdiagm(0 => Sy_p * H' * uᵧ)
        Kz = spdiagm(0 => Sz_p * H' * uᵧ)
        diagW = diag(blockdiag(W[1], W[2], W[3]))
        new_diagW = [val != 0 ? 1.0 / val : 1.0 for val in diagW]
        Wꜝ = spdiagm(0 => new_diagW)
        return ConvectionOps{N}((Cx, Cy, Cz), (Kx, Ky, Kz), G, H, Wꜝ, V, size)
    end
end


"""
    struct AdvectionVecOps{N} <: AbstractOperators where N

Struct representing a collection of advection operators for a vector field.
"""
struct AdvectionVecOps{N} <: AbstractOperators
    C :: NTuple{N, SparseMatrixCSC{Float64, Int}}
    K :: NTuple{N, SparseMatrixCSC{Float64, Int}}
    V :: NTuple{N, SparseMatrixCSC{Float64, Int}}
    H :: NTuple{N, SparseMatrixCSC{Float64, Int}}
    size :: NTuple{N, Int}
end

# 1D Case
function AdvectionVecOps(cap::NTuple{1},size,uₒ,uᵧ)
    nx = size[1] + 2

    Dx_m = ẟ_m(nx)
    Dx_p = δ_p(nx)
    Sx_m = Σ_m(nx)
    Sx_p = Σ_p(nx)

    H = cap[1].A[1]*Dx_m - Dx_m*cap[1].B[1]

    Cx = Dx_p * spdiagm(0 => (Sx_m * cap[1].A[1] * uₒ[1])) * Sx_m

    Kx = spdiagm(0 => Sx_p * H' * uᵧ[1])

    return AdvectionVecOps{1}((Cx,), (Kx,), (cap[1].V,), (H,), size)
end

# 2D Case
function AdvectionVecOps(cap::NTuple{2}, size, uₒ, uᵧ)
    nxu, nyu = size[1]+2, size[2]+1
    nxv, nyv = size[1]+1, size[2]+2

    # For u-velocity (faces in x-direction)
    Dx_mu, Dy_mu = kron(I(nyu), ẟ_m(nxu)), kron(ẟ_m(nyu), I(nxu)) # Dxx-; Dxy-
    Dx_pu, Dy_pu = kron(I(nyu), δ_p(nxu)), kron(δ_p(nyu), I(nxu)) # Dxx+; Dxy+
    Sx_mu, Sy_mu = kron(I(nyu), Σ_m(nxu)), kron(Σ_m(nyu), I(nxu)) # Sxx-; Sxy-
    Sx_pu, Sy_pu = kron(I(nyu), Σ_p(nxu)), kron(Σ_p(nyu), I(nxu)) # Sxx+; Sxy+

    # For v-velocity (faces in y-direction)
    Dx_mv, Dy_mv = kron(I(nyv), ẟ_m(nxv)), kron(ẟ_m(nyv), I(nxv)) # Dyx-; Dyy-
    Dx_pv, Dy_pv = kron(I(nyv), δ_p(nxv)), kron(δ_p(nyv), I(nxv)) # Dyx+; Dyy+
    Sx_mv, Sy_mv = kron(I(nyv), Σ_m(nxv)), kron(Σ_m(nyv), I(nxv)) # Syx-; Syy-
    Sx_pv, Sy_pv = kron(I(nyv), Σ_p(nxv)), kron(Σ_p(nyv), I(nxv)) # Syx+; Syy+
    
    Hu = [cap[1].A[1]*Dx_mu - Dx_mu*cap[1].B[1] ; cap[1].A[2]*Dy_mu - Dy_mu*cap[1].B[2]]
    Hv = [cap[2].A[1]*Dx_mv - Dx_mv*cap[2].B[1] ; cap[2].A[2]*Dy_mv - Dy_mv*cap[2].B[2]]

    Cx1 = Dx_pu * spdiagm(0 => (Sx_mu * cap[1].A[1] * uₒ[1])) * Sx_mu
    Cx2 = Dy_pu * spdiagm(0 => (Sx_mu * cap[1].A[2] * uₒ[1])) * Sy_mu
    Cy1 = Dx_pv * spdiagm(0 => (Sy_mv * cap[2].A[1] * uₒ[2])) * Sx_mv
    Cy2 = Dy_pv * spdiagm(0 => (Sy_mv * cap[2].A[2] * uₒ[2])) * Sy_mv

    Cx = Cx1 + Cx2
    Cy = Cy1 + Cy2

    Kx = spdiagm(0 => Sx_pu * Hu' * vcat(uᵧ[1], uᵧ[1]))
    Ky = spdiagm(0 => Sy_pv * Hv' * vcat(uᵧ[2], uᵧ[2]))

    return AdvectionVecOps{2}((Cx, Cy), (Kx, Ky), (cap[1].V, cap[2].V), (Hu, Hv), size)
end


struct SpaceTimeOps{N} <: AbstractOperators where N
    G :: SparseMatrixCSC{Float64, Int}
    H :: SparseMatrixCSC{Float64, Int}
    Wꜝ :: SparseMatrixCSC{Float64, Int}
    V :: SparseMatrixCSC{Float64, Int}
    size :: NTuple{N, Int}
end

function SpaceTimeOps(A, B, V, W, size)
    N = length(size)
    if N == 1
        nx = size[1]
        G = ẟ_m(nx) * B[1] # Gérer le periodicity
        H = A[1]*ẟ_m(nx) - ẟ_m(nx)*B[1]
        diagW = diag(blockdiag(W[1]))
        new_diagW = [val != 0 ? 1.0 / val : 1.0 for val in diagW]
        Wꜝ = spdiagm(0 => new_diagW)
    elseif N == 2
        nx, ny = size[1], size[2]
        Dx_m = kron(I(ny), ẟ_m(nx))
        Dy_m = kron(ẟ_m(ny), I(nx))
        G = [Dx_m * B[1]; Dy_m * B[2]]
        H = [A[1]*Dx_m - Dx_m*B[1]; A[2]*Dy_m - Dy_m*B[2]]
        diagW = diag(blockdiag(W[1], W[2]))
        new_diagW = [val != 0 ? 1.0 / val : 1.0 for val in diagW]
        Wꜝ = spdiagm(0 => new_diagW)
    elseif N == 3
        nx, ny, nz = size[1], size[2], size[3]
        Dx_m = kron(I(nz), kron(I(ny), ẟ_m(nx)))
        Dy_m = kron(I(nz), kron(ẟ_m(ny), I(nx)))
        Dz_m = kron(ẟ_m(nz), kron(I(ny), I(nx)))
        G = [Dx_m * B[1]; Dy_m * B[2]; Dz_m * B[3]]
        H = [A[1]*Dx_m - Dx_m*B[1]; A[2]*Dy_m - Dy_m*B[2]; A[3]*Dz_m - Dz_m*B[3]]
        diagW = diag(blockdiag(W[1], W[2], W[3]))
        new_diagW = [val != 0 ? 1.0 / val : 1.0 for val in diagW]
        Wꜝ = spdiagm(0 => new_diagW)
    end

    return SpaceTimeOps{N}(G, H, Wꜝ, V, size)
end


# Volume redefinition
function volume_redefinition!(capacity::Capacity{1}, operator::DiffusionOps)
    pₒ = [capacity.C_ω[i][1] for i in 1:length(capacity.C_ω)]
    pᵧ = [capacity.C_γ[i][1] for i in 1:length(capacity.C_ω)]
    p = vcat(pₒ, pᵧ)
    grad = ∇(operator, p)
    W_new = [grad[i] * capacity.W[1][i,i] for i in 1:length(grad)]
    W_new = spdiagm(0 => W_new)

    pₒ = [(capacity.C_ω[i][1]^2)/2 for i in 1:length(capacity.C_ω)]
    pᵧ = [(capacity.C_γ[i][1]^2)/2 for i in 1:length(capacity.C_ω)]

    p = vcat(pₒ, pᵧ)
    grad = ∇(operator, p)

    qω = vcat(grad)
    qγ = vcat(grad)

    div = ∇_(operator, qω, qγ)
    V_new = spdiagm(0 => div)

    capacity.W = (W_new,)
    capacity.V = V_new
end