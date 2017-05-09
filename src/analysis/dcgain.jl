dcgain(s::StateSpace{Val{:siso}}) = _dcgain(s)[1]
dcgain(s::StateSpace{Val{:mimo}}) = _dcgain(s)
dcgain(s::LtiSystem)              = dcgain(ss(s))

_dcgain{S}(s::StateSpace{S,Val{:cont}}) = s.D - s.C/s.A*s.B
_dcgain{S}(s::StateSpace{S,Val{:disc}}) = s.C/(I - s.A)*s.B + s.D
