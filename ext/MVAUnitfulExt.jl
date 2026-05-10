module MVAUnitfulExt

import MinimumVarianceAnalysis: mva_eigen, _field, _quantity_mva_eigen
using Unitful: Quantity, ustrip, unit, dimension
using Unitful: V, m, T

_field_dimension(q) = _field(dimension(q), dimension(V / m), dimension(T))

function mva_eigen(B::AbstractMatrix{Q}; field = nothing, kwargs...) where {Q <: Quantity}
    return _quantity_mva_eigen(B, unit(Q), ustrip, _field_dimension; field, kwargs...)
end

end
