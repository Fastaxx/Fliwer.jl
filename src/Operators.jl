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
        Cx = spdiagm(ẟ_m(nx) * Σ_p(nx) * A[1] * uₒ) * Σ_p(nx)
        H = A[1]*ẟ_m(nx) - ẟ_m(nx)*B[1]
        Kx = Σ_m(nx) * spdiagm(0 => H' * uᵧ)
        return ConvectionOps{N}((Cx,), (Kx,), size)
    elseif N == 2
        nx, ny = size[1], size[2]
        Dx_m = kron(I(ny), ẟ_m(nx))
        Dy_m = kron(ẟ_m(ny), I(nx))
        Cx = spdiagm(0 => Dx_m * Σ_p(nx) * A[1] * uₒ[1]) * Σ_p(nx) + spdiagm(0 => Dy_m * Σ_p(nx) * A[2] * uₒ[2]) * Σ_p(ny)
        Cy = spdiagm(0 => Dx_m * Σ_p(ny) * A[1] * uₒ[1]) * Σ_p(nx) + spdiagm(0 => Dy_m * Σ_p(ny) * A[2] * uₒ[2]) * Σ_p(ny)
        H = [A[1]*Dx_m - Dx_m*B[1]; A[2]*Dy_m - Dy_m*B[2]]
        Kx = Σ_m(nx) * spdiagm(0 => H' * uᵧ)
        Ky = Σ_m(ny) * spdiagm(0 => H' * uᵧ)
        return ConvectionOps{N}((Cx, Cy), (Kx, Ky), size)
    elseif N == 3
        nx, ny, nz = size[1], size[2], size[3]
        Dx_m = kron(I(nz), kron(I(ny), ẟ_m(nx)))
        Dy_m = kron(I(nz), kron(ẟ_m(ny), I(nx)))
        Dz_m = kron(ẟ_m(nz), kron(I(ny), I(nx)))
        Cx = spdiagm(0 => Dx_m * Σ_p(nx) * A[1] * uₒ[1]) * Σ_p(nx) + spdiagm(0 => Dy_m * Σ_p(nx) * A[2] * uₒ[2]) * Σ_p(ny) + spdiagm(0 => Dz_m * Σ_p(nx) * A[3] * uₒ[3]) * Σ_p(nz)
        Cy = spdiagm(0 => Dx_m * Σ_p(ny) * A[1] * uₒ[1]) * Σ_p(nx) + spdiagm(0 => Dy_m * Σ_p(ny) * A[2] * uₒ[2]) * Σ_p(ny) + spdiagm(0 => Dz_m * Σ_p(ny) * A[3] * uₒ[3]) * Σ_p(nz)
        Cz = spdiagm(0 => Dx_m * Σ_p(nz) * A[1] * uₒ[1]) * Σ_p(nx) + spdiagm(0 => Dy_m * Σ_p(nz) * A[2] * uₒ[2]) * Σ_p(ny) + spdiagm(0 => Dz_m * Σ_p(nz) * A[3] * uₒ[3]) * Σ_p(nz)
        H = [A[1]*Dx_m - Dx_m*B[1]; A[2]*Dy_m - Dy_m*B[2]; A[3]*Dz_m - Dz_m*B[3]]
        Kx = Σ_m(nx) * spdiagm(0 => H' * uᵧ)
        Ky = Σ_m(ny) * spdiagm(0 => H' * uᵧ)
        Kz = Σ_m(nz) * spdiagm(0 => H' * uᵧ)
        return ConvectionOps{N}((Cx, Cy, Cz), (Kx, Ky, Kz), size)
    end
end