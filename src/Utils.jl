@inline CI(a...) = CartesianIndex(a...) # abréviation pour la création d'objets CartesianIndex

"""
    CIj(j, I, jj)
Replace jᵗʰ component of CartesianIndex with k
"""
CIj(j, I::CartesianIndex{d}, k) where d = CI(ntuple(i -> i == j ? k : I[i], d))

"""
    δ(i, N::Int)
    δ(i, I::CartesianIndex{N}) where {N}

    Return a CartesianIndex of dimension `N` which is one at index `i` and zero elsewhere.
"""
δ(i, ::Val{N}) where N = CI(ntuple(j -> j == i ? 1 : 0, N))
δ(i, I::CartesianIndex{N}) where N = δ(i, Val{N}())


"""
    inside(a)

    Return CartesianIndices range excluding a single layer of cells on all boundaries.
"""
@inline inside(a::AbstractArray; buff=1) = CartesianIndices(map(ax -> first(ax) + buff : last(ax) - buff, axes(a)))


"""
    @inside <expr>

Simple macro to automate efficient loops over cells excluding ghosts. For example,

    @inside p[I] = sum(loc(0,I))

becomes

    @loop p[I] = sum(loc(0,I)) over I ∈ inside(p)

See [`@loop`](@ref).
"""
macro inside(ex)
    # Make sure it's a single assignment
    @assert ex.head == :(=) && ex.args[1].head == :(ref)
    a,I = ex.args[1].args[1:2]
    return quote # loop over the size of the reference
        WaterLily.@loop $ex over $I ∈ inside($a)
    end |> esc
end


"""
    inside_u(dims,j)

Return CartesianIndices range excluding the ghost-cells on the boundaries of
a _vector_ array on face `j` with size `dims`.
"""
function inside_u(dims::NTuple{N},j) where {N}
    CartesianIndices(ntuple( i-> i==j ? (3:dims[i]-1) : (2:dims[i]), N))
end
@inline inside_u(dims::NTuple{N}) where N = CartesianIndices((map(i->(2:i-1),dims)...,1:N))
@inline inside_u(u::AbstractArray) = CartesianIndices(map(i->(2:i-1),size(u)[1:end-1]))
splitn(n) = Base.front(n),last(n)
size_u(u) = splitn(size(u))

"""
    loc(i,I) = loc(Ii)

Location in space of the cell at CartesianIndex `I` at face `i`.
Using `i=0` returns the cell center s.t. `loc = I`.
"""
@inline loc(i,I::CartesianIndex{N},T=Float32) where N = SVector{N,T}(I.I .- 1.5 .- 0.5 .* δ(i,I).I)
@inline loc(Ii::CartesianIndex,T=Float32) = loc(last(Ii),Base.front(Ii),T)
Base.last(I::CartesianIndex) = last(I.I)
Base.front(I::CartesianIndex) = CI(Base.front(I.I))