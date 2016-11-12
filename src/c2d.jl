module Discretization

abstract Method

immutable ZOH <: Method
end

@compat function (m::ZOH){T}(s::StateSpace{T,Siso{true}}, Ts::Real)
  A, B, C, D  = s.A, s.B, s.C, s.D
  ny, nu      = size(s)
  nx          = numstates(s)
  M     = expm([A*Ts  B*Ts;
          zeros(nu, nx + nu)])
  Ad    = M[1:nx, 1:nx]
  Bd    = M[1:nx, nx+1:nx+nu]
  Cd    = C
  Dd    = D
  x0map = [speye(nx) spzeros(nx, nu)]
  Ad, Bd, Cd, Dd, Ts, x0map
end

immutable FOH <: Method
end

@compat function (m::FOH){T}(s::StateSpace{T,Continuous{true}}, Ts::Real)
  A, B, C, D  = s.A, s.B, s.C, s.D
  ny, nu      = size(s)
  nx          = numstates(s)
  M     = expm([A*Ts B*Ts zeros(nx, nu);
          zeros(nu, nx + nu) eye(nu);
          zeros(nu, nx + 2*nu)])
  M1    = M[1:nx, nx+1:nx+nu]
  M2    = M[1:nx, nx+nu+1:nx+2*nu]
  Ad    = M[1:nx, 1:nx]
  Bd    = Ad*M2 + M1 - M2
  Cd    = C
  Dd    = D + C*M2
  x0map = [eye(nx) -M2]
  Ad, Bd, Cd, Dd, Ts, x0map
end

immutable GeneralizedBilinear{T} <: Method
  α::T
  @compat function (::Type{GeneralizedBilinear}){T}(α::T = 0.5)
    # TODO: Any assertions needed?
    new{T}(α)
  end
end

@compat function (m::GeneralizedBilinear){T}(s::StateSpace{T,Siso{true}}, Ts::Real)
  A, B, C, D  = s.A, s.B, s.C, s.D
  ny, nu      = size(s)
  nx          = numstates(s)
  α           = m.α
  ima   = eye(nx) - α*Ts*A
  Ad    = ima\(eye(nx) + (1.0-α)*Ts*A)
  Bd    = ima\(Ts*B)
  Cd    = (ima.'\C.').'
  Dd    = D + α*(C*Bd)
  x0map = [eye(nx) zeros(nx, nu)]
  Ad, Bd, Cd, Dd, Ts, x0map
end

immutable ForwardEuler <: Method
end

@compat function (m::ForwardEuler){T}(s::StateSpace{T,Siso{true}}, Ts::Real)
  A, B, C, D  = s.A, s.B, s.C, s.D
  ny, nu      = size(s)
  nx          = numstates(s)
  Ad    = eye(nx) + Ts*A
  Bd    = Ts*B
  Cd    = C
  Dd    = D
  x0map = [eye(nx) zeros(nx, nu)]
  Ad, Bd, Cd, Dd, Ts, x0map
end

immutable BackwardEuler <: Method
end

@compat function (m::BackwardEuler){T}(s::StateSpace{T,Siso{true}}, Ts::Real)
  A, B, C, D  = s.A, s.B, s.C, s.D
  ny, nu      = size(s)
  nx          = numstates(s)
  ima   = eye(nx) - Ts*A
  Ad    = ima\eye(nx)
  Bd    = ima\(Ts*B)
  Cd    = (ima.'\C.').'
  Dd    = D + C*Bd
  x0map = [eye(nx) zeros(nx, nu)]
  Ad, Bd, Cd, Dd, Ts, x0map
end

# TODO: Discretization.Method implementations for s::LtiSystem
# TODO: A brainstorming regarding support for any callable objects:
#       - Do we need to specialize c2d on RationalTF and ZeroPoleGain, as well,
#         and document the method, accordingly, or,
#       - Do we just require the method to use s::LtiSystem, no matter what, but
#         use a type assertion inside c2d such that the return value of the method
#         is a SISO/MIMO discrete system (whichever is suitable for the function
#         call)?

end

"""
    c2d(s, Ts[, method])

Convert the continuous system `s` into a discrete system with sample time
`Ts`, using the provided discretization method with zero-order-hold as default.

`method` can be any callable object with a signature `method(s, Ts)`.

`c2d` also supports `do ... end` block calls in the form

    c2d(s, Ts) do s, Ts
      # Your transformation based on `s` and `Ts`
    end

For a system in state space form, returns the discretized system as well as a
matrix `x0map` that transforms the initial conditions to the discrete domain by
`x0_discrete = x0map*[x0; u0]`.

For a system in transfer function form, only the discretized transfer function
is returned.

Discretization methods:

- zero-order-hold                 (`Discretization.ZOH()`)
- first-order-hold                (`Discretization.FOH()`)
- bilinear transform              (`Discretization.Bilinear()`)
- forward Euler                   (`Discretization.ForwardEuler()`)
- backward Euler                  (`Discretization.BackwardEuler()`)
- generalized bilinear transform  (`Discretization.Bilinear(α::Real)`)

The generalized bilinear transform uses the parameter α and is based on [1].

- [1] G. Zhang, X. Chen, and T. Chen, Digital redesign via the generalized
bilinear transformation, Int. J. Control, vol. 82, no. 4, pp. 741-754, 2009.
"""
function c2d(s::StateSpace{Siso{true},Continuous{true}}, Ts::Real,
  method = Discretization.ZOH())
  @assert Ts > zero(Ts) && !isinf(Ts) "c2d: Ts must be a positive number"
  Ad, Bd, Cd, Dd, Ts, x0map = method(s, Ts)
  ss(Ad, Bd, Cd, Dd[1], Ts), x0map
end

function c2d(s::StateSpace{Siso{false},Continuous{true}}, Ts::Real,
  method = Discretization.ZOH())
  @assert Ts > zero(Ts) && !isinf(Ts) "c2d: Ts must be a positive number"
  Ad, Bd, Cd, Dd, Ts, x0map = method(s, Ts)
  ss(Ad, Bd, Cd, Dd, Ts), x0map
end

c2d{T}(s::LtiSystem{T,Continuous{true}}, Ts::Real,
  method = Discretization.ZOH())
  @assert Ts > zero(Ts) && !isinf(Ts) "c2d: Ts must be a positive number"
  method(s, Ts)::LtiSystem{T,Continuous{false}}
end

c2d{T}(method::Function, s::StateSpace{T,Continuous{true}}, Ts::Real) = c2d(s, Ts, method)
c2d{T}(method::Function, s::LtiSystem{T,Continuous{true}}, Ts::Real) = c2d(s, Ts, method)

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
