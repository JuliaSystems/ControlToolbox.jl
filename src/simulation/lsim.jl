# Discrete State Space

"""
    `lsim(sys, u, t[, x0, method])`

Calculate the time response of system `sys` to input `u`. If `x0` is ommitted,
a zero vector is used.

Continuous time systems are discretized before simulation. By default, the
discretization method is chosen based on the smoothness of the input signal. Optionally, the
`method` parameter can be specified as either `:zoh` or `:foh`.
"""
function lsim{T1<:Real,T2<:AbstractVector}(s::DSs, u::Array{T1}, t::T2=[],
    x0=zeros(T1, s.nx, 1))
  @assert length(x0) == s.nx string("length(x0) must equal number of states of sys")
  @assert size(u,2) == s.nu string("Number of input signals need to match system")
  N, nu = size(u)
  x = Array(T1, s.nx, N)
  for i=1:N
      x[:,i] = x0
      x0 = s.A * x0 + s.B * u[i,:].'
  end
  y = s.C*x + s.D*(u.')
  return y.', x.'
end

function lsim{T1<:AbstractVector}(s::DSs, u::Function, t::T1=[], x0::AbstractVector=zeros(Int8, s.nx, 1))
  @assert length(x0) == s.nx string("length(x0) must equal number of states of sys")
  x = Array(Float64, s.nx, N)
  U = Array(Float64, N, 1)
  for i=1:N
    x[:,i] = x0
    U[i] = u(x0)
    x0 = s.A * x0 + s.B * U[i]
  end
  y = s.C*x + s.D*(U.')
  return y', x', t, U
end

# Discrete SISO Transfer Function
function lsim{T1<:Real,T2<:AbstractVector}(s::ControlCore.DSisoTf, u::Array{T1},
    t::T2=[], x0...)

  # zeropad b and a if necessary
  b = numvec(s)
  a = denvec(s)
  lengthb = length(b)
  lengtha = length(a)
  order = max(lengthb,lengtha)
  b = copy!(zeros(eltype(b),order),order-lengthb+1,b,1)
  a = copy!(zeros(eltype(a),order),order-lengtha+1,a,1)

  length(x0) > 0 ? filt(b, a, u, x0[1]) :
                   filt(b, a, u)
end

#  Discrete MIMO Transfer Function
function lsim{T1<:Real,T2<:AbstractVector}(s::ControlCore.DMimo, u::Array{T1},
    t::T2=[], x0...)
  @assert s.nu == size(u,2) string("Number of input signals need to match system")
  U = Array(Array{T1,1},s.ny,s.nu)
  for i = 1:s.ny
    for j = 1:s.nu
      U[i,j] = u[:,j]
    end
  end
  length(x0) > 0 ?
    map((x,u,x0_)->lsim(x,u,t,x0_...), ControlCore.getmatrix(s), U, x0[1])  :
    map((x,u)->lsim(x,u,t), ControlCore.getmatrix(s), U)
end

# Continuous Versions
function lsim{T1<:Real,T2<:AbstractVector}(s::CSs, u::Array{T1}, t::T2=[],
    x0=zeros(T1, s.nx, 1), method::Symbol=:zoh)
  @assert length(t) > 1 string("length(t) must be at least 2")
  Ts = t[2]-t[1]
  ds, x0map = c2d(s, Ts, method)
  x0 = x0map*[x0; u[1,:].']
  lsim(ds, u, t, x0)
end

function lsim{T1<:Real,T2<:AbstractVector}(s::ControlCore.CSisoTf, u::Array{T1},
    t::T2, method::Symbol=:zoh)
  @assert length(t) > 1 string("length(t) must be at least 2")
  Ts = t[2]-t[1]
  lsim(c2d(s, Ts, method), u, t)
end

function lsim{T1<:Real,T2<:AbstractVector}(s::ControlCore.CMimo, u::Array{T1},
    t::T2=[], method::Symbol=:zoh)
  @assert size(u,2) == s.nu string("Number of input signals need to match system")
  @assert length(t) > 1 string("length(t) must be at least 2")
  Ts = t[2]-t[1]
  lsim(c2d(s, Ts, method), u, t)
end

function lsim{T1<:AbstractVector,T2<:AbstractVector}(s::CSs, u::Function, t::T1,
    x0::T2=zeros(Int8,s.nx, 1), method::Symbol=:zoh)
  @assert length(x0) == s.nx string("length(x0) must equal number of states of sys")
  @assert length(t) > 2
  Ts = t[2]-t[1]
  ds, x0map = c2d(sys, Ts, method)
  x0 = x0map*[x0; u(1,x0)]
  lsim(ds, u, t, x0)
end
