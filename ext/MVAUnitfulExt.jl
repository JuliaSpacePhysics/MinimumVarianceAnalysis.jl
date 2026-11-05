module MVAUnitfulExt

import MinimumVarianceAnalysis: mva_eigen, _quantity_mva_eigen
using Unitful: Quantity, ustrip, unit, dimension, V, m

mva_eigen(B::AbstractMatrix{<:Quantity}; kwargs...) =
    _quantity_mva_eigen(B, Base.Fix1(ustrip, unit(first(B))), dimension(first(B)) == dimension(V / m); kwargs...)

end
