module MVAUnitfulExt

import MinimumVarianceAnalysis: mva_eigen, _quantity_mva_eigen
using Unitful: Quantity, ustrip, dimension, V, m

mva_eigen(B::AbstractMatrix{<:Quantity}; kwargs...) =
    _quantity_mva_eigen(B, ustrip(B), dimension(first(B)) == dimension(V / m); kwargs...)

end
