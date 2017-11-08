dcgain(sys::LTISystems.LtiSystem) = map(real, sys(ω = 0.))
