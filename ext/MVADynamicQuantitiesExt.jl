module MVADynamicQuantitiesExt

import MinimumVarianceAnalysis: mva_eigen, _field, _quantity_mva_eigen
using DynamicQuantities: Dimensions, Quantity, UnionAbstractQuantity, ustrip, dimension

const ELECTRIC_DIMENSION = Dimensions(length = 1, mass = 1, time = -3, current = -1)
const MAGNETIC_DIMENSION = Dimensions(mass = 1, time = -2, current = -1)

_field_dimension(q) = _field(dimension(q), ELECTRIC_DIMENSION, MAGNETIC_DIMENSION)

function _scale(B)
    return Quantity(1.0, dimension(B))
end

function mva_eigen(B::AbstractMatrix{<:UnionAbstractQuantity}; field = nothing, kwargs...)
    return _quantity_mva_eigen(B, _scale(B), ustrip, _field_dimension; field, kwargs...)
end

end
