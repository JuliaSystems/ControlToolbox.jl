"""
    d2c(s)

Recovers the continuous time which was sampled with a zero-order-hold to obtain
the discrete time system `s`.
"""
function d2c{S}(s::StateSpace{S,Val{:disc}})
  A, B, C, D = s.A, s.B, s.C, s.D
  for λ in eig(A)[1]
    println(λ)
    if λ < zero(λ)
      error("d2c cannot be applied to systems with negative eigenvalues")
    end
  end
#  @assert ~any(eig(A) .< zero(eltype(A))) "d2c cannot be applied to systems with negative eigenvalues"
  ny, nu = size(s)
  nx = s.nx
  M = logm([A  B; zeros(nu, nx) eye(nu)]) / s.Ts
  Ac = M[1:nx, 1:nx]
  Bc = M[1:nx, nx+1:nx+nu]
  Cc = C
  Dc = D
  ss(Ac, Bc, Cc, Dc)
end

d2c{S}(s::LtiSystem{S,Val{:disc}}) = d2c(ss(s))
