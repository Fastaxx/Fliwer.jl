abstract type Operators end

# Définition de la structure DiffusionOps
struct DiffusionOps <: Operators
    G::SparseMatrixCSC{Float64, Int}
    H::SparseMatrixCSC{Float64, Int}
    Wꜝ::SparseMatrixCSC{Float64, Int}
end

# Opérateurs élémentaires
function ẟ_m(n::Int, periodicity::Bool=false) D = spdiagm(0 => ones(n), -1 => ones(n-1)); D[n, n] = 0.0; if periodicity; D[1, n-1] = -1.0; D[n, 1] = 1.0; end; D end
function δ_p(n::Int, periodicity::Bool=false) D = spdiagm(0 => -ones(n), 1 => ones(n-1)); D[n, n] = 0.0; if periodicity; D[1, n-1] = -1.0; D[n, 1] = 1.0; end; D end
function Σ_m(n::Int, periodicity::Bool=false) D = 0.5 * spdiagm(0 => ones(n), -1 => ones(n-1)); D[n, n] = 0.0; if periodicity; D[1, n-1] = 0.5; D[n, 1] = 0.5; end; D end
function Σ_p(n::Int, periodicity::Bool=false) D = 0.5 * spdiagm(0 => ones(n), 1 => ones(n-1)); D[n, n] = 0.0; if periodicity; D[1, n-1] = 0.5; D[n, 1] = 0.5; end; D end
function I(n::Int) spdiagm(0 => ones(n)) end

# Fonction pour créer DiffusionOps à partir de A, B, V, W
function DiffusionOps(A, B, V, W)
    N = length(A) 
    if N == 1
        nx = size(A[1], 1)
        G = ẟ_m(nx) * B[1] # Gérer le periodicity
        H = A[1]*ẟ_m(nx) - ẟ_m(nx)*B[1]
        Wꜝ = spdiagm(0 => diagm(W) .!= 0 ? 1 ./ diagm(W) : 1.0)
    elseif N == 2
        nx, ny = size(A[1], 1), size(A[2], 1)
        Dx_m = kron(I(ny), ẟ_m(nx))
        Dy_m = kron(ẟ_m(ny), I(nx))
        G = [Dx_m * B[1]; Dy_m * B[2]]
        H = [A[1]*Dx_m - Dx_m*B[1]; A[2]*Dy_m - Dy_m*B[2]]
        Wd = spdiagm(0 => blockdiag(W[1], W[2]) .!= 0 ? 1 ./ blockdiag(W[1], W[2]) : 1.0) 
    elseif N == 3
        nx, ny, nz = size(A[1], 1), size(A[2], 1), size(A[3], 1)
        Dx_m = kron(I(nz), kron(I(ny), ẟ_m(nx)))
        Dy_m = kron(I(nz), kron(ẟ_m(ny), I(nx)))
        Dz_m = kron(ẟ_m(nz), kron(I(ny), I(nx)))
        G = [Dx_m * B[1]; Dy_m * B[2]; Dz_m * B[3]]
        H = [A[1]*Dx_m - Dx_m*B[1]; A[2]*Dy_m - Dy_m*B[2]; A[3]*Dz_m - Dz_m*B[3]]
        Wd = spdiagm(0 => blockdiag(W[1], W[2], W[3]) .!= 0 ? 1 ./ blockdiag(W[1], W[2], W[3]) : 1.0)  
    end

    return DiffusionOps(G, H, Wd)
end

struct ConvectionOps <: Operators
    Conv::SparseMatrixCSC{Float64, Int}
end