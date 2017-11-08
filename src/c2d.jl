module Discretization

using Compat
using LTISystems
using LTISystems: LtiSystem, StateSpace, TransferFunction

# Method abstraction
@compat abstract type Method end

# Zero-Order-Hold
immutable ZOH <: Method
end

@compat function (m::ZOH)(A::AbstractMatrix, B::AbstractMatrix, C::AbstractMatrix,
  D::AbstractMatrix, Ts::Real)
  ny, nu      = size(D)
  nx          = size(A,1)
  M     = expm([A*Ts  B*Ts;
          zeros(nu, nx + nu)])
  Ad    = M[1:nx, 1:nx]
  Bd    = M[1:nx, nx+1:nx+nu]
  Cd    = C
  Dd    = D
  x0map = [speye(nx) spzeros(nx, nu)]
  Ad, Bd, Cd, Dd, Ts, x0map
end
@compat function (m::ZOH)(s::StateSpace{Val{:siso},Val{:cont}}, Ts::Real)
  Ad, Bd, Cd, Dd, Ts, x0map = m(s.A, s.B, s.C, s.D, Ts)
  ss(Ad, Bd, Cd, Dd[1], Ts), x0map
end
@compat function (m::ZOH)(s::StateSpace{Val{:mimo},Val{:cont}}, Ts::Real)
  Ad, Bd, Cd, Dd, Ts, x0map = m(s.A, s.B, s.C, s.D, Ts)
  ss(Ad, Bd, Cd, Dd, Ts), x0map
end
@compat (m::ZOH){T}(s::TransferFunction{Val{T},Val{:cont}}, Ts::Real) =
  tf(m(ss(s), Ts)[1])
@compat (m::ZOH){T}(s::LtiSystem{Val{T},Val{:cont}}, Ts::Real)  =
  m(ss(s), Ts)[1]

# First-Order-Hold
immutable FOH <: Method
end

@compat function (m::FOH)(A::AbstractMatrix, B::AbstractMatrix, C::AbstractMatrix,
  D::AbstractMatrix, Ts::Real)
  ny, nu      = size(D)
  nx          = size(A,1)
  M     = expm([A*Ts B*Ts zeros(nx, nu);
          zeros(nu, nx + nu) eye(nu);
          zeros(nu, nx + 2*nu)])
  M1    = M[1:nx, nx+1:nx+nu]
  M2    = M[1:nx, nx+nu+1:nx+2*nu]
  Ad    = M[1:nx, 1:nx]
  Bd    = Ad*M2 + M1 - M2
  Cd    = C
  Dd    = D + C*M2
  x0map = [speye(nx) -M2]
  Ad, Bd, Cd, Dd, Ts, x0map
end
@compat function (m::FOH)(s::StateSpace{Val{:siso},Val{:cont}}, Ts::Real)
  Ad, Bd, Cd, Dd, Ts, x0map = m(s.A, s.B, s.C, s.D, Ts)
  ss(Ad, Bd, Cd, Dd[1], Ts), x0map
end
@compat function (m::FOH)(s::StateSpace{Val{:mimo},Val{:cont}}, Ts::Real)
  Ad, Bd, Cd, Dd, Ts, x0map = m(s.A, s.B, s.C, s.D, Ts)
  ss(Ad, Bd, Cd, Dd, Ts), x0map
end
@compat (m::FOH){T}(s::TransferFunction{Val{T},Val{:cont}}, Ts::Real) =
  tf(m(ss(s), Ts)[1])
@compat (m::FOH){T}(s::LtiSystem{Val{T},Val{:cont}}, Ts::Real)  =
  m(ss(s), Ts)[1]

# Generalized Bilinear Transformation
immutable Bilinear{T<:Real} <: Method
  α::T
  @compat function (::Type{Bilinear}){T}(α::T = 0.5)
    @assert α ≥ 0. && α ≤ 1. "Bilinear: α must be between 0 and 1"
    new{T}(α)
  end
end

@compat function (m::Bilinear)(A::AbstractMatrix, B::AbstractMatrix,
  C::AbstractMatrix, D::AbstractMatrix, Ts::Real)
  ny, nu      = size(D)
  nx          = size(A,1)
  α           = m.α
  ima   = I - α*Ts*A
  Ad    = ima\(I + (1.0-α)*Ts*A)
  Bd    = ima\(Ts*B)
  Cd    = (ima.'\C.').'
  Dd    = D + α*(C*Bd)
  x0map = [speye(nx) spzeros(nx, nu)]
  Ad, Bd, Cd, Dd, Ts, x0map
end
@compat function (m::Bilinear)(s::StateSpace{Val{:siso},Val{:cont}}, Ts::Real)
  Ad, Bd, Cd, Dd, Ts, x0map = m(s.A, s.B, s.C, s.D, Ts)
  ss(Ad, Bd, Cd, Dd[1], Ts), x0map
end
@compat function (m::Bilinear)(s::StateSpace{Val{:mimo},Val{:cont}}, Ts::Real)
  Ad, Bd, Cd, Dd, Ts, x0map = m(s.A, s.B, s.C, s.D, Ts)
  ss(Ad, Bd, Cd, Dd, Ts), x0map
end
@compat (m::Bilinear){T}(s::TransferFunction{Val{T},Val{:cont}}, Ts::Real)  =
  tf(m(ss(s), Ts)[1])
@compat (m::Bilinear){T}(s::LtiSystem{Val{T},Val{:cont}}, Ts::Real)   =
  m(ss(s), Ts)[1]

# Forward Euler
immutable ForwardEuler <: Method
end

@compat function (m::ForwardEuler)(A::AbstractMatrix, B::AbstractMatrix,
  C::AbstractMatrix, D::AbstractMatrix, Ts::Real)
  ny, nu      = size(D)
  nx          = size(A,1)
  Ad    = I + Ts*A
  Bd    = Ts*B
  Cd    = C
  Dd    = D
  x0map = [speye(nx) spzeros(nx, nu)]
  Ad, Bd, Cd, Dd, Ts, x0map
end
@compat function (m::ForwardEuler)(s::StateSpace{Val{:siso},Val{:cont}}, Ts::Real)
  Ad, Bd, Cd, Dd, Ts, x0map = m(s.A, s.B, s.C, s.D, Ts)
  ss(Ad, Bd, Cd, Dd[1], Ts), x0map
end
@compat function (m::ForwardEuler)(s::StateSpace{Val{:mimo},Val{:cont}}, Ts::Real)
  Ad, Bd, Cd, Dd, Ts, x0map = m(s.A, s.B, s.C, s.D, Ts)
  ss(Ad, Bd, Cd, Dd, Ts), x0map
end
@compat (m::ForwardEuler){T}(s::TransferFunction{Val{T},Val{:cont}}, Ts::Real)  =
  tf(m(ss(s), Ts)[1])
@compat (m::ForwardEuler){T}(s::LtiSystem{Val{T},Val{:cont}}, Ts::Real)   =
  m(ss(s), Ts)[1]

# Backward Euler
immutable BackwardEuler <: Method
end

@compat function (m::BackwardEuler)(A::AbstractMatrix, B::AbstractMatrix,
  C::AbstractMatrix, D::AbstractMatrix, Ts::Real)
  ny, nu      = size(D)
  nx          = size(A,1)
  ima   = I - Ts*A
  Ad    = ima\eye(nx)
  Bd    = ima\(Ts*B)
  Cd    = (ima.'\C.').'
  Dd    = D + C*Bd
  x0map = [speye(nx) spzeros(nx, nu)]
  Ad, Bd, Cd, Dd, Ts, x0map
end
@compat function (m::BackwardEuler)(s::StateSpace{Val{:siso},Val{:cont}}, Ts::Real)
  Ad, Bd, Cd, Dd, Ts, x0map = m(s, Ts)
  ss(Ad, Bd, Cd, Dd[1], Ts), x0map
end
@compat function (m::BackwardEuler)(s::StateSpace{Val{:mimo},Val{:cont}}, Ts::Real)
  Ad, Bd, Cd, Dd, Ts, x0map = m(s, Ts)
  ss(Ad, Bd, Cd, Dd, Ts), x0map
end
@compat (m::BackwardEuler){T}(s::TransferFunction{Val{T},Val{:cont}}, Ts::Real)  =
  tf(m(ss(s), Ts)[1])
@compat (m::BackwardEuler){T}(s::LtiSystem{Val{T},Val{:cont}}, Ts::Real)   =
  m(ss(s), Ts)[1]

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
function c2d{T}(s::StateSpace{T,Val{:cont}}, Ts::Real,
  method = Discretization.ZOH())
  @assert Ts > zero(Ts) && !isinf(Ts) "c2d: Ts must be a positive number"
  sys, x0map = method(s, Ts)
  return sys::StateSpace{T,Val{:disc}}, x0map::AbstractMatrix
end
function c2d{T}(s::LtiSystem{T,Val{:cont}}, Ts::Real,
  method = Discretization.ZOH())
  @assert Ts > zero(Ts) && !isinf(Ts) "c2d: Ts must be a positive number"
  method(s, Ts)::LtiSystem{T,Val{:disc}}
end
c2d{T}(method::Function, s::LtiSystem{T,Val{:cont}}, Ts::Real) = c2d(s, Ts, method)
