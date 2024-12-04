"""
    abstract type AbstractOperators

An abstract type representing a collection of operators.
"""
abstract type AbstractOperators end

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
        Cx = δ_p(nx) * spdiagm(0 => diag(Σ_m(nx) * A[1] * uₒ[1])) * Σ_m(nx)
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
        Cx = Dx_p * spdiagm(0 => diag(Sx_m * A[1] * uₒ[1])) * Sx_m
        Cy = Dy_p * spdiagm(0 => diag(Sy_m * A[2] * uₒ[2])) * Sy_m
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
        Cx = Dx_p * spdiagm(0 => diag(Sx_m * A[1] * uₒ[1])) * Sx_m
        Cy = Dy_p * spdiagm(0 => diag(Sy_m * A[2] * uₒ[2])) * Sy_m
        Cz = Dz_p * spdiagm(0 => diag(Sz_m * A[3] * uₒ[3])) * Sz_m
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