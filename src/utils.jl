cross3(x, y) = cross(SV3(x), SV3(y))

is_right_handed(v1, v2, v3) = cross3(v1, v2) ⋅ v3 > 0

@views function is_right_handed(F::Eigen)
    vs = F.vectors
    v1 = vs[:, 1]
    v2 = vs[:, 2]
    v3 = vs[:, 3]
    return is_right_handed(v1, v2, v3)
end

_dimnum(x, query = nothing) = 1

"""
    transform(A, mat::AbstractMatrix; dim=nothing, query=nothing)

Transform `A` into a new coordinate system using transformation matrix `mat` along the `dim` dimension (time).
"""
function transform(A, mat::AbstractMatrix; dim = nothing, query = nothing)
    dim = @something dim _dimnum(A, query)
    return dim == 1 ? A * mat : mat' * A
end

transform(A, mat::Eigen; kw...) = transform(A, mat.vectors; kw...)

# Deprecated alias for backwards compatibility
rotate(A, mat; kw...) = transform(A, mat; kw...)

"""
    normal(F::Eigen; field=:B)

Return the boundary normal eigenvector from an MVA result.

- `field=:B` (MVAB): minimum variance direction (last eigenvector)
- `field=:E` (MVAE): maximum variance direction (first eigenvector)
"""
@views normal(F::Eigen; field = :B) = field == :E ? F.vectors[:, 1] : F.vectors[:, end]
