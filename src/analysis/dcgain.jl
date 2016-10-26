dcgain(s::LtiSystem) = dcgain(ss(s))

function _dcgain{S}(s::StateSpace{S,Continuous{true}})
  return -s.C/s.A*s.B + s.D
end

function _dcgain{S}(s::StateSpace{S,Continuous{false}})
  return s.C/(speye(size(s.A)) - s.A)*s.B + s.D
end

dcgain(s::StateSpace{Siso{true}}) = _dcgain(s)[1]
dcgain(s::StateSpace{Siso{false}}) = _dcgain(s)
