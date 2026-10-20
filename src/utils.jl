cross3(x, y) = cross(SV3(x), SV3(y))

is_right_handed(v1, v2, v3) = cross3(v1, v2) ⋅ v3 > 0

@views function is_right_handed(F::Eigen)
    vs = F.vectors
    v1 = vs[:, 1]
    v2 = vs[:, 2]
    v3 = vs[:, 3]
    return is_right_handed(v1, v2, v3)
end

"""
    normal(F::Eigen; field=:B)

Return the boundary normal eigenvector from an MVA result.

- `field=:B` (MVAB): minimum variance direction (last eigenvector)
- `field=:E` (MVAE): maximum variance direction (first eigenvector)
"""
@views normal(F::Eigen; field = :B) = field == :E ? F.vectors[:, 1] : F.vectors[:, end]
