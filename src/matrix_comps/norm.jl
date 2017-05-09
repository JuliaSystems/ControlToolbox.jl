"""
`norm(sys[, Val{p}])`

Compute the `p`-norm of the system `sys`. `p` can be either `2` or `Inf`
(default is 2)
"""
function norm{S}(sys::StateSpace{S,Val{:cont}}, ::Type{Val{2}})
  A, B, C, D = sys.A, sys.B, sys.C, sys.D
  if !isstable(sys) || any(D .!= 0)
    return Inf
  else
    Q = clyap(A, B*B')
    return sqrt(trace(C*Q*C'))
  end
end
function norm{S}(sys::StateSpace{S,Val{:disc}}, ::Type{Val{2}})
  A, B, C, D = sys.A, sys.B, sys.C, sys.D
  if !isstable(sys)
    return Inf
  else
    Q = dlyap(A, B*B')
    return sqrt(trace(C*Q*C' + D*D'))
  end
end
norm(sys::LtiSystem) = norm(sys,Val{2})
norm(sys::LtiSystem,n::Int) = norm(sys,Val{n})
norm{T}(sys::LtiSystem,::Type{Val{T}}) = norm(ss(sys),Val{T})
