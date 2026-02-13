using Unitful: Quantity, ustrip, unit, dimension
using Unitful: V, m, T

function _field(unit)
    d = dimension(unit)
    d == dimension(V / m) && return :E
    d == dimension(T) && return :B
    return nothing
end

function mva_eigen(B::AbstractMatrix{Q}; field = nothing, kwargs...) where {Q <: Quantity}
    field = @something field _field(unit(Q)) :B
    F = mva_eigen(ustrip(B); field, kwargs...)
    return Eigen(F.values * unit(Q)^2, F.vectors)
end