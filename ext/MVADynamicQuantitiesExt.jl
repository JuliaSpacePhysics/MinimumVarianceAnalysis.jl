module MVADynamicQuantitiesExt

import MinimumVarianceAnalysis: mva_eigen, _quantity_mva_eigen
using DynamicQuantities: Dimensions, UnionAbstractQuantity, ustrip, dimension

const ELECTRIC_DIMENSION = Dimensions(length = 1, mass = 1, time = -3, current = -1)

mva_eigen(B::AbstractMatrix{<:UnionAbstractQuantity}; kwargs...) =
    _quantity_mva_eigen(B, ustrip.(B), dimension(first(B)) == ELECTRIC_DIMENSION; kwargs...)

end
