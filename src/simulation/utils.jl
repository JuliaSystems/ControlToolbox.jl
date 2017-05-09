function _default_time_vector(s::LtiSystem, Tf::Real=-1)
  Ts = _default_Ts(s)
  if Tf < zero(Int)
    Tf = 99*Ts
  end
  return 0:Ts:Tf
end

function _default_Ts(s::LtiSystem)
  if isdiscrete(s)
    Ts = s.Ts
  elseif !isstable(s)
    Ts = 0.05
  else
    ps = pole(s)
    r = minimum([abs(real(ps));0])
  if r == zero(Float64)
    r = 1.0
  end
    Ts = 0.07/r
  end
  return Ts
end
