"""
`covar(sys, W)`

Calculate the stationary covariance of an LTI model `sys`, driven by Gaussian
white noise of covariance `W`
"""
function covar{S}(sys::StateSpace{S,Val{:cont}}, W::StridedMatrix)
  A, B, C, D = sys.A, sys.B, sys.C, sys.D
  if size(B,2) != size(W, 1) || size(W, 1) != size(W, 2)
    error("W must be a square matrix the same size as `sys.B` columns")
  end
  if !isstable(sys) || any(D .!= 0)
    return Inf
  else
    Q = clyap(A, B*W*B')
    return C*Q*C'
  end
end
function covar{S}(sys::StateSpace{S,Val{:disc}}, W::StridedMatrix)
  A, B, C, D = sys.A, sys.B, sys.C, sys.D
  if size(B,2) != size(W, 1) || size(W, 1) != size(W, 2)
    error("W must be a square matrix the same size as `sys.B` columns")
  end
  if !isstable(sys)
    return Inf
  else
    Q = dlyap(A, B*W*B')
    return C*Q*C' + D*W*D'
  end
end
covar(sys::LtiSystem,W::StridedMatrix) = covar(ss(sys),W)
