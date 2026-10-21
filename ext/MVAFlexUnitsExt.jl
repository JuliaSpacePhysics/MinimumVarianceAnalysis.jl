module MVAFlexUnitsExt

import MinimumVarianceAnalysis: mva_eigen, _quantity_mva_eigen
using FlexUnits: Dimensions, QuantUnion, ustrip, dimension

# m kg s A K cd mol; positional because the field names changed in FlexUnits 0.6
const ELECTRIC_DIMENSION = Dimensions(1, 1, -3, -1, 0, 0, 0)

mva_eigen(B::AbstractMatrix{<:QuantUnion}; kwargs...) =
    _quantity_mva_eigen(B, ustrip.(B), dimension(first(B)) == ELECTRIC_DIMENSION; kwargs...)

end
