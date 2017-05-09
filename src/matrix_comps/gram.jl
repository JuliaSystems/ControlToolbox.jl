"""
`gram(sys, Val{opt})`

Compute the grammian of state-space system `sys`. If `opt` is `:c`, computes the
controllability grammian. If `opt` is `:o`, computes the observability
grammian.
"""
function gram{S,T}(sys::StateSpace{S,Val{T}}, ::Type{Val{:c}})
  if !isstable(sys)
    error("gram only valid for stable systems")
  else
    return lyap(sys.A, sys.B*sys.B',Val{T})
  end
end
function gram{S,T}(sys::StateSpace{S,Val{T}}, ::Type{Val{:o}})
  if !isstable(sys)
    error("gram only valid for stable systems")
  else
    return lyap(sys.A', sys.C'*sys.C,Val{T})
  end
end
function gram{T}(sys::StateSpace, ::Type{Val{T}})
  error("opt must be either :c for controllability grammian, or :o for observability grammian")
end
gram(sys::StateSpace, T::Symbol) = gram(sys,Val{T})
