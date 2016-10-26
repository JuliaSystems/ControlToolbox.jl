
isstable{S}(s::LtiSystem{S,Continuous{false}}) = maximum(abs(poles(s))) < 1

isstable{S}(s::LtiSystem{S,Continuous{true}}) = maximum(real(poles(s))) < 0

isstable(s::LtiSystem{Siso{false}}) = map(isstable, getmatrix(s))

# TODO: How to handle pole-zero cancellation

# Legacy
#isstable{T<:ControlCore.DSiso}(s::T) = maximum(abs([poles(s); 0])) < 1
#isstable{T<:ControlCore.CSiso}(s::T) = maximum(real([poles(s); 0])) < 0
#isstable{T<:ControlCore.MimoSystem}(s::T) = map(isstable, getmatrix(s))
