"""
    c2d(s, Ts, method[, α])

Convert the continuous system `s` into a discrete system with sample time
`Ts`, using the provided discretization method with zero-order-hold as default.

For a system in state space form, returns the discretized system as well as a
matrix `x0map` that transforms the initial conditions to the discrete domain by
`x0_discrete = x0map*[x0; u0]`.

For a system in transfer function form, only the discretized transfer function
is returned.

Discretization methods:

- zero-order-hold     (`:zoh`)
- first-order-hold    (`:foh`)
- bilinear transform (`:bilinear` or `:tustin`)
- forward euler (`:euler`  or `:forward_diff`)
- backward euler (`:backward_diff`)
- generalized bilinear transform (`:gbt`)

The generalized bilinear transform uses the additional parameter α and is based on [1].

- [1] G. Zhang, X. Chen, and T. Chen, Digital redesign via the generalized
bilinear transformation, Int. J. Control, vol. 82, no. 4, pp. 741-754, 2009.
"""
function c2d{S}(s::LtiSystem{S,Continuous{true}}, Ts::Real, method::Symbol=:zoh,
                α::Real=zero(Float64))
  if method == :zoh
    return c2dzoh(s, Ts)
  elseif method == :foh
    return c2dfoh(s, Ts)
  elseif method == :bilinear || method == :tustin
    return c2dgbt(s, Ts, 0.5)
  elseif method == :euler || method == :forward_diff
    return c2dforward(s, Ts)
  elseif method == :backward_diff
    return c2dbackward(s, Ts)
  elseif method == :gbt
    return c2dgbt(s, Ts, α)
  else
    error("Unsupported method: ", method)
  end
end

# Internal methods
function c2dzoh{S}(s::StateSpace{S,Continuous{true}}, Ts::Real)
  A, B, C, D = s.A, s.B, s.C, s.D
  ny, nu = size(s)
  nx = s.nx
  M = expm([A*Ts  B*Ts;
            zeros(nu, nx + nu)])
  Ad = M[1:nx, 1:nx]
  Bd = M[1:nx, nx+1:nx+nu]
  Cd = C
  Dd = D
  x0map = [speye(nx) spzeros(nx, nu)]
  ss(Ad, Bd, Cd, Dd, Ts), x0map
end

c2dzoh{S}(s::LtiSystem{S,Continuous{true}}, Ts::Real) = c2dzoh(ss(s),Ts)[1]

function c2dfoh{S}(s::StateSpace{S,Continuous{true}}, Ts::Real)
  A, B, C, D = s.A, s.B, s.C, s.D
  ny, nu = size(s)
  nx = s.nx
  M = expm([A*Ts B*Ts zeros(nx, nu);
       zeros(nu, nx + nu) eye(nu);
       zeros(nu, nx + 2*nu)])
  M1 = M[1:nx, nx+1:nx+nu]
  M2 = M[1:nx, nx+nu+1:nx+2*nu]
  Ad = M[1:nx, 1:nx]
  Bd = Ad*M2 + M1 - M2
  Cd = C
  Dd = D + C*M2
  x0map = [eye(nx) -M2]
  ss(Ad, Bd, Cd, Dd, Ts), x0map
end

c2dfoh{S}(s::LtiSystem{S,Continuous{true}}, Ts::Real) = c2dfoh(s,Ts)[1]

function c2dgbt{S}(s::StateSpace{S,Continuous{true}}, Ts::Real,
                                      α::Real=zero(Float64))
  A, B, C, D = s.A, s.B, s.C, s.D
  ny, nu = size(s)
  nx = s.nx
  ima = eye(nx) - α*Ts*A
  Ad = ima\(eye(nx) + (1.0-α)*Ts*A)
  Bd = ima\(Ts*B)
  Cd = (ima.'\C.').'
  Dd = D + α*(C*Bd)
  x0map = [eye(nx) zeros(nx, nu)]
  ss(Ad, Bd, Cd, Dd, Ts), x0map
end

c2dgbt{S}(s::LtiSystem{S,Continuous{true}}, Ts::Real, α::Real=zero(Float64)) = c2dgbt(s,Ts,α)[1]

function c2dforward{S}(s::StateSpace{S,Continuous{true}}, Ts::Real)
  A, B, C, D = s.A, s.B, s.C, s.D
  ny, nu = size(s)
  nx = s.nx
  Ad = eye(nx) + Ts*A
  Bd = Ts*B
  Cd = C
  Dd = D
  x0map = [eye(nx) zeros(nx, nu)]
  ss(Ad, Bd, Cd, Dd, Ts), x0map
end

c2dforward{S}(s::LtiSystem{S,Continuous{true}}, Ts::Real) = c2dforward(ss(s),Ts)[1]

function c2dbackward{S}(s::StateSpace{S,Continuous{true}}, Ts::Real)
  A, B, C, D = s.A, s.B, s.C, s.D
  ny, nu = size(s)
  nx = s.nx
  ima = eye(nx) - Ts*A
  Ad = ima\eye(nx)
  Bd = ima\(Ts*B)
  Cd = (ima.'\C.').'
  Dd = D + C*Bd
  x0map = [eye(nx) zeros(nx, nu)]
  ss(Ad, Bd, Cd, Dd, Ts), x0map
end

c2dbackward{S}(s::LtiSystem{S,Continuous{true}}, Ts::Real) = c2dbackward(ss(s),Ts)[1]

# Zhai, Guisheng, et al. "An extension of generalized bilinear transformation for digital redesign."
# International Journal of Innovative Computing, Information and Control 8.6 (2012): 4071-4081.
# http://www.ijicic.org/icmic10-ijicic-03.pdf



# Legacy
# DSs = Union{ControlCore.DSisoSs,ControlCore.DMimoSs}
# CSs = Union{ControlCore.CSisoSs,ControlCore.CMimoSs}

# # catch all
# function c2d{T1<:ControlCore.SisoSystem,T2<:Real,T3<:Real}(
#     s::T1, Ts::T2, method::Symbol=:zoh, α::T3=zero(Float64))
#   #tf(c2d(tf(s), Ts, method, α)[1])
# end
#
# function c2d{T1<:ControlCore.SisoSystem,T2<:Real,T3<:Real}(
#     s::ControlCore.MimoSystem{T1}, Ts::T2,
#     method::Symbol=:zoh, α::T3=zero(Float64))
#   #zpk(c2d(zpk2ss(s), Ts, method, α)[1])
# end
#
#
# function c2d{T1<:ControlCore.CSisoRational,M1,T2<:Real,T3<:Real}(
#     s::ControlCore.CMimo{T1,M1}, Ts::T2,
#     method::Symbol=:zoh, α::T3=zero(Float64))
#   #tf(c2d(ss(s), Ts, method, α)[1])
# end
#
# function c2d{T1<:ControlCore.CSisoRational,T2<:Real,T3<:Real}(s::T1, Ts::T2,
#     method::Symbol=:zoh, α::T3=zero(Float64))
#   #tf(c2d(ss(s), Ts, method, α)[1])
# end
#
# function c2d{T1<:ControlCore.CSisoZpk,M1,T2<:Real,T3<:Real}(
#     s::ControlCore.CMimo{T1,M1}, Ts::T2,
#     method::Symbol=:zoh, α::T3=zero(Float64))
#   #zpk(c2d(ss(s), Ts, method, α)[1])
# end
#
# function c2d{T1<:ControlCore.CSisoZpk,T2<:Real,T3<:Real}(s::T1, Ts::T2,
#     method::Symbol=:zoh, α::T3=zero(Float64))
#   #zpk(c2d(ss(s), Ts, method, α)[1])
# end

# function c2d{T1<:Real,T2<:Real}(s::CSs, Ts::T1, method::Symbol=:zoh, α::T2=zero(Float64))
#   if method == :zoh
#     return c2dzoh(s, Ts)
#   elseif method == :foh
#     return c2dfoh(s, Ts)
#   elseif method == :bilinear || method == :tustin
#     return c2dgbt(s, Ts, 0.5)
#   elseif method == :euler || method == :forward_diff
#     return c2dforward(s, Ts)
#   elseif method == :backward_diff
#     return c2dbackward(s, Ts)
#   elseif method == :gbt
#     return c2dgbt(s, Ts, α)
#   else
#     error("Unsupported method: ", method)
#   end
# end

# Internal methods

# function c2dzoh{T1<:Real}(s::CSs, Ts::T1)
#   A, B, C, D = s.A, s.B, s.C, s.D
#   ny, nu = size(s)
#   nx = s.nx
#   M = expm([A*Ts  B*Ts;
#         zeros(nu, nx + nu)])
#   Ad = M[1:nx, 1:nx]
#   Bd = M[1:nx, nx+1:nx+nu]
#   Cd = C
#   Dd = D
#   x0map = [eye(nx) zeros(nx, nu)]
#   ss(Ad, Bd, Cd, Dd, Ts), x0map
# end
#
# function c2dfoh{T1<:Real}(s::CSs, Ts::T1)
#   A, B, C, D = s.A, s.B, s.C, s.D
#   ny, nu = size(s)
#   nx = s.nx
#   M = expm([A*Ts B*Ts zeros(nx, nu);
#        zeros(nu, nx + nu) eye(nu);
#        zeros(nu, nx + 2*nu)])
#   M1 = M[1:nx, nx+1:nx+nu]
#   M2 = M[1:nx, nx+nu+1:nx+2*nu]
#   Ad = M[1:nx, 1:nx]
#   Bd = Ad*M2 + M1 - M2
#   Cd = C
#   Dd = D + C*M2
#   x0map = [eye(nx) -M2]
#   ss(Ad, Bd, Cd, Dd, Ts), x0map
# end
#
# function c2dgbt{T1<:Real,T2<:Real}(s::CSs, Ts::T1, α::T2=zero(Float64))
#   A, B, C, D = s.A, s.B, s.C, s.D
#   ny, nu = size(s)
#   nx = s.nx
#   ima = eye(nx) - α*Ts*A
#   Ad = ima\(eye(nx) + (1.0-α)*Ts*A)
#   Bd = ima\(Ts*B)
#   Cd = (ima.'\C.').'
#   Dd = D + α*(C*Bd)
#   x0map = [eye(nx) zeros(nx, nu)]
#   ss(Ad, Bd, Cd, Dd, Ts), x0map
# end
#
# function c2dforward{T1<:Real}(s::CSs, Ts::T1)
#   A, B, C, D = s.A, s.B, s.C, s.D
#   ny, nu = size(s)
#   nx = s.nx
#   Ad = eye(nx) + Ts*A
#   Bd = Ts*B
#   Cd = C
#   Dd = D
#   x0map = [eye(nx) zeros(nx, nu)]
#   ss(Ad, Bd, Cd, Dd, Ts), x0map
# end
#
# function c2dbackward{T1<:Real}(s::CSs, Ts::T1)
#   A, B, C, D = s.A, s.B, s.C, s.D
#   ny, nu = size(s)
#   nx = s.nx
#   ima = eye(nx) - Ts*A
#   Ad = ima\eye(nx)
#   Bd = ima\(Ts*B)
#   Cd = (ima.'\C.').'
#   Dd = D + C*Bd
#   x0map = [eye(nx) zeros(nx, nu)]
#   ss(Ad, Bd, Cd, Dd, Ts), x0map
# end
