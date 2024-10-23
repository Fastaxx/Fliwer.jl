abstract type AbstractOperators end

# Définition de la structure DiffusionOps
struct DiffusionOps{N} <: AbstractOperators where N
    G::SparseMatrixCSC{Float64, Int}
    H::SparseMatrixCSC{Float64, Int}
    Wꜝ::SparseMatrixCSC{Float64, Int}
    V::SparseMatrixCSC{Float64, Int}
    size::NTuple{N, Int}
end

# Opérateurs élémentaires
function ẟ_m(n::Int, periodicity::Bool=false) D = spdiagm(0 => ones(n), -1 => -ones(n-1)); D[n, n] = 0.0; if periodicity; D[1, n-1] = -1.0; D[n, 1] = 1.0; end; D end
function δ_p(n::Int, periodicity::Bool=false) D = spdiagm(0 => -ones(n), 1 => ones(n-1)); D[n, n] = 0.0; if periodicity; D[1, n-1] = -1.0; D[n, 1] = 1.0; end; D end
function Σ_m(n::Int, periodicity::Bool=false) D = 0.5 * spdiagm(0 => ones(n), -1 => ones(n-1)); D[n, n] = 0.0; if periodicity; D[1, n-1] = 0.5; D[n, 1] = 0.5; end; D end
function Σ_p(n::Int, periodicity::Bool=false) D = 0.5 * spdiagm(0 => ones(n), 1 => ones(n-1)); D[n, n] = 0.0; if periodicity; D[1, n-1] = 0.5; D[n, 1] = 0.5; end; D end
function I(n::Int) spdiagm(0 => ones(n)) end

# Fonction pour créer DiffusionOps à partir de A, B, V, W
function DiffusionOps(A, B, V, W, size)
    N = length(size)
    if N == 1
        nx = size[1]
        G = ẟ_m(nx) * B[1] # Gérer le periodicity
        H = A[1]*ẟ_m(nx) - ẟ_m(nx)*B[1]
        Wꜝ = spdiagm(0 => float.([val !=0 ? 1/val : 1 for val in diag(W)]))
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

struct ConvectionOps <: AbstractOperators
    Conv::SparseMatrixCSC{Float64, Int}
end