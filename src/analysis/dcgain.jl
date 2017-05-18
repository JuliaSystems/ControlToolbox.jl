dcgain(sys::SystemsBase.LtiSystem) = map(real, sys(ω = 0.))
