dcgain(s::StateSpace{Siso{true}}) = _dcgain(s)[1]
dcgain(s::StateSpace{Siso{false}}) = _dcgain(s)
dcgain(s::LtiSystem) = dcgain(ss(s))

_dcgain{S}(s::StateSpace{S,Continuous{true}}) = s.D - s.C/s.A*s.B
_dcgain{S}(s::StateSpace{S,Continuous{false}}) = s.C/(speye(size(s.A)) - s.A)*s.B + s.D
