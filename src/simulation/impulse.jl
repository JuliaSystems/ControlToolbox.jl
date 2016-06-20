"""

    `impulse(s)`

Calculate the impulse response of system `s`.


If the final time `Tf` or time vector `t` is not provided, one is calculated based on the system pole
locations.

`y` has size `(length(t), ny, nu)`, `x` has size `(length(t), nx, nu)`
"""
# Discrete State Space
function impulse{T1<:AbstractVector}(s::DSs, t::T1)
  nu = numinputs(s)
  ny = numoutputs(s)
  nx = numstates(s)
  u = zeros(length(t), nu)
  x = zeros(length(t), nx, nu)
  y = zeros(length(t), ny, nu)
  x0s = zeros(nx, nu)
  for i=1:nu
    u[1,:] = zeros(1, nu)
    u[1,i] = 1/samplingtime(s)
    y[:,:,i], x[:,:,i]  = lsim(s, u, t, x0s)
  end
  return y, t, x
end

# Discrete SISO Transfer Function
function impulse{T1<:AbstractVector}(s::ControlCore.DSisoTf, t::T1)
  u = zeros(length(t), 1)
  u[1] = 1/s.Ts
  y = lsim(s, u, t)
  return y, t
end

#  Discrete MIMO Transfer Function
function impulse{T1<:AbstractVector}(s::ControlCore.DMimo, t::T1)
  nu = numinputs(s)
  ny = numoutputs(s)
  u = zeros(length(t), nu)
  y = zeros(length(t), ny, nu)
  for i=1:nu
    u[1,:] = zeros(1,nu)
    u[1,i] = 1/samplingtime(s)
    y[:,:,i] = lsim(s, u, t)
  end
  return y, t
end

# Continuous Versions
function impulse{T1<:AbstractVector}(s::CSs, t::T1)
  # impulse response equivalent to unforced response of
  # ss(A, 0, C, 0) with x0 = B.
  nu = numinputs(s)
  ny = numoutputs(s)
  nx = numstates(s)
  u = zeros(length(t), nu)
  x = zeros(length(t), nx, nu)
  y = zeros(length(t), ny, nu)
  imp_s = ss(s.A, zeros(nx, 1), s.C, zeros(ny, 1))
  for i=1:nu
    y[:,:,i], x[:,:,i]  = lsim(imp_s, u, t, s.B[:,i], :zoh)
  end
  return y, t, x
end

function impulse{T1<:AbstractVector}(s::ControlCore.CSisoTf, t::T1)
  impulse(ss(s), t)[1:2]
end

function impulse{T1<:AbstractVector}(s::ControlCore.CMimo, t::T1)
  impulse(ss(s), t)[1:2]
end

impulse(s::ControlCore.LtiSystem, Tf::Real) = impulse(s, _default_time_vector(s, Tf))
impulse(s::ControlCore.LtiSystem) = impulse(s, _default_time_vector(s))
