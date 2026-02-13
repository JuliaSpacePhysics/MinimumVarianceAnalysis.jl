module MVADimensionalDataExt
import MinimumVarianceAnalysis: mva_eigen, _dimnum
using DimensionalData: AbstractDimArray, dimnum, hasdim, TimeDim

function _dimnum(x::AbstractDimArray, query = nothing)
    # If query is nothing, choose the time dimension or the first dimension
    return if isnothing(query)
        hasdim(x, TimeDim) ? dimnum(x, TimeDim) : 1
    else
        dimnum(x, query)
    end
end

function mva_eigen(A::AbstractDimArray; dim = nothing, query = nothing, kw...)
    dim = @something dim _dimnum(A, query)
    return mva_eigen(parent(A); dim, kw...)
end
end
