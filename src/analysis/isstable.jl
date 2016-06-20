
isstable{T<:ControlCore.DSiso}(s::T) = maximum(abs([poles(s); 0])) < 1

isstable{T<:ControlCore.CSiso}(s::T) = maximum(real([poles(s); 0])) < 0

isstable{T<:ControlCore.MimoSystem}(s::T) = map(isstable, getmatrix(s))
# TODO: How to handle pole-zero cancellation
