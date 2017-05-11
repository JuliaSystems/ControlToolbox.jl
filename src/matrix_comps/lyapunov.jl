"""
`clyap(A, Q)`

Compute the solution `X` to the discrete Lyapunov equation
`AX + XA' + Q = 0`.
"""
# TODO: Change code by SLICOT version
function clyap{T<:BlasFloat}(A::StridedMatrix{T}, Q::StridedMatrix{T})
  lhs = kron(speye(size(A)...), A) + kron(conj(A),speye(size(A)...))
  x = -lhs\reshape(Q, prod(size(Q)), 1)
  return reshape(x, size(Q))
end
clyap{T1<:Real,T2<:Real}(A::StridedMatrix{T1}, Q::StridedMatrix{T2}) =
  clyap(float(A), float(Q))


"""
`dlyap(A, Q)`

Compute the solution `X` to the discrete Lyapunov equation
`AXA' - X + Q = 0`.
"""
# TODO: Change code by SLICOT version
function dlyap{T<:BlasFloat}(A::StridedMatrix{T}, Q::StridedMatrix{T})
  lhs = kron(conj(A), A)
  lhs = speye(size(lhs)...) - lhs
  x = lhs\reshape(Q, prod(size(Q)), 1)
  return reshape(x, size(Q))
end
dlyap{T1<:Real,T2<:Real}(A::StridedMatrix{T1}, Q::StridedMatrix{T2}) =
  dlyap(float(A), float(Q))

lyap{T<:BlasFloat}(A::StridedMatrix{T}, Q::StridedMatrix{T},::Type{Val{:cont}}) =
  clyap(A,Q)

lyap{T<:BlasFloat}(A::StridedMatrix{T}, Q::StridedMatrix{T},::Type{Val{:disc}}) =
  dlyap(A,Q)
