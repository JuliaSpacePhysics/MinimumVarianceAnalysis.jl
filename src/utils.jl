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

function rotate(A, mat::AbstractMatrix; dim = nothing, query = nothing)
    dim = @something dim _dimnum(A, query)
    return dim == 1 ? A * mat : mat' * A
end

rotate(A, mat::Eigen; kw...) = rotate(A, mat.vectors; kw...)