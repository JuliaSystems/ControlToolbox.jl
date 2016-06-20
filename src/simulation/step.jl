"""
    `step(s[,t])`

Calculate the step response of system `s`.

If the final time `t` or time vector `t` is not provided, one is calculated based on the system pole
locations.
"""


step{T1<:AbstractVector}(s::ControlCore.DSisoTf, t::T1) = lsim(s, ones(length(t)), t), t
step{T1<:AbstractVector}(s::ControlCore.CSisoTf, t::T1) = lsim(s, ones(length(t)), t, :zoh), t

function step{T1<:AbstractVector}(s::ControlCore.DMimo, t::T1)
  nu = numinputs(s)
  ny = numoutputs(s)
  u  = zeros(length(t), nu)
  y  = zeros(length(t), ny, nu)
  for i=1:nu
    u        = zeros(length(t), nu)
    u[:,i]   =  ones(length(t))
    y[:,:,i] = lsim(s, u, t)
  end
  return y, t
end

function step{T1<:AbstractVector}(s::ControlCore.CMimo, t::T1)
  nu = numinputs(s)
  ny = numoutputs(s)
  u  = zeros(length(t), nu)
  y  = zeros(length(t), ny, nu)
  for i=1:nu
    u        = zeros(length(t), nu)
    u[:,i]   = ones(length(t))
    y[:,:,i] = lsim(s, u, t, :zoh)
  end
  return y, t
end

# for state-space return also state information
function step{T1<:AbstractVector}(s::Union{DSs,CSs}, t::T1)
  nu = numinputs(s)
  ny = numoutputs(s)
  nx = numstates(s)
  u  = zeros(length(t), nu)
  x  = zeros(length(t), nx, nu)
  y  = zeros(length(t), ny, nu)
  for i=1:nu
    u                  = zeros(length(t), nu)
    u[:,i]             =  ones(length(t))
    y[:,:,i], x[:,:,i] = lsim(s, u, t, zeros(nx, 1), :zoh)
  end
  return y, t, x
end

step(s::ControlCore.LtiSystem, Tf::Real) = step(s, default_time_vector(s, Tf))
step(s::ControlCore.LtiSystem) = step(s, default_time_vector(s))
