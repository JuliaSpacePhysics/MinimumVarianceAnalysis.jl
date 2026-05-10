module MVADimensionalDataExt
import MinimumVarianceAnalysis: _dimnum
using DimensionalData: AbstractDimArray, dimnum, TimeDim

_dimnum(x::AbstractDimArray, dim = nothing) = dimnum(x, @something dim TimeDim)
end
