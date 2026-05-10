module MVAFlexUnitsExt

import MinimumVarianceAnalysis: mva_eigen, _field, _quantity_mva_eigen
using FlexUnits: Dimensions, QuantUnion, dimension, quantity, ustrip

const ELECTRIC_DIMENSION = Dimensions(length = 1, mass = 1, time = -3, current = -1)
const MAGNETIC_DIMENSION = Dimensions(mass = 1, time = -2, current = -1)

_raw_values(B) = ustrip.(B)

_field_dimension(q) = _field(dimension(q), ELECTRIC_DIMENSION, MAGNETIC_DIMENSION)

function _scale(B)
    return quantity(1.0, dimension(first(B)))
end

function mva_eigen(B::AbstractMatrix{<:QuantUnion}; field = nothing, kwargs...)
    return _quantity_mva_eigen(B, _scale(B), _raw_values, _field_dimension; field, kwargs...)
end

end
