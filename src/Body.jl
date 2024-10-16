
"""
    AbstractBody

Embedded body Abstract type. Any `AbstractBody` subtype must implement

d = sdf(body::AbstractBody, x, t=0)

and

d,n,V = measure(body::AbstractBody, x, t=0, fastd²=Inf)

where `d` is the signed distance from `x` to the body at time `t`,
and `n` & `V` are the normal and velocity vectors implied at `x`.
A fast-approximate method can return `≈d,zero(x),zero(x)` if `d^2>fastd²`.
"""
abstract type AbstractBody end

"""
    NoBody

Use for a simulation without any body.
"""
struct NoBody <: AbstractBody end

"""
    measure_sdf!(a::AbstractArray, body::AbstractBody, t::Number = 0; kwargs...)

Remplit l'array `a` avec le champ de distance signée (SDF) par rapport à `body` au temps `t`.

# Arguments
- `a` : Array multidimensionnel à remplir avec les valeurs SDF.
- `body` : Géométrie immergée utilisée pour calculer le SDF.
- `t` : Temps actuel de la simulation (optionnel, par défaut `0`).
- `kwargs` : Arguments supplémentaires passés à la fonction `sdf`.
"""
function measure_sdf!(a::AbstractArray, body::AbstractBody, t::Number = 0; kwargs...)
    # Définir les plages d'indices internes pour chaque dimension, excluant une couche de bordure
    internal_ranges = [first(ax) + 1 : last(ax) - 1 for ax in axes(a)]
    
    # Générer les indices cartésiens pour les cellules internes
    internal_indices = CartesianIndices(internal_ranges...)
    
    # Itérer sur les indices internes et remplir l'array avec les valeurs SDF
    for I in internal_indices
        # Calculer les coordonnées physiques à partir des indices
        # Supposons une grille avec un espacement unitaire et une origine à zéro
        coords = Tuple((I[d] - 1) for d in 1:ndims(a))
        
        # Calculer la distance signée à l'emplacement `coords` pour le temps `t`
        a[I] = sdf(body, coords, t; kwargs...)
    end
end

