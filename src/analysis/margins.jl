"""
`pm, wpm = phasemargin(sys; deg=true/false)`

Compute the phase margin and corresponding unity gain crossover
frequency for LTI system `sys` or data from bode-plot `mag`, `phase`, `w`.
By definition, unity gain feedback with phase `pm` would make the closed
loop system unstable. By default `pm` is given in radians.
"""
function phasemargin(sys::TransferFunction{Val{:siso},Val{:cont}}; deg::Bool=false)
  # For sys = 1/(s+1)^2*(s+10)^2*(s+100)^4/100^4/10 the output from phasemargin is off by a noticable ammount. Haven't been able to replicate or find the problem, might be because of instabilities in roots() though...

  Z, P = sys.num[1], sys.den[1]

  # To find w180, first solve |H(jw)| == 1
  rZ, iZ, rP, iP = polysplit(Z)..., polysplit(P)...
  wctmp = roots(rZ^2 + iZ^2 - rP^2 - iP^2)

  # wc must be real and positive
  wc = real(wctmp[(imag(wctmp) .== 0) & (real(wctmp) .>= 0)])
  isempty(wc) && (return Inf, NaN)

  # Evaluate at the valid frequencies and pick the one where the phase is closest to an odd multiple of π
  pm = angle(vec(freqresp(sys, wc)[1])) # note: -π < pm < π
  ipm = findmax(abs(pm))[2]

  return (pm[ipm]<0 ? π + pm[ipm] : pm[ipm] - π)*(deg ? 360/2π : 1), wc[ipm]
end

# For other representations, convert to transfer function
phasemargin(sys::LtiSystem{Val{:siso}}; deg::Bool=false) = phasemargin(tf(sys), deg=deg)


"""
`gm, wgm = gainmargin(sys; dB=true/false)`

Compute the gain margin and corresponding 180 degree phase
crossing for the LTI system `sys`. The gain margin is the smallest
ammount of (pure gain) negative feedback that would make the
closed loop system unstable. By default the result is given
in linear units (not dB).
"""
function gainmargin(sys::TransferFunction{Val{:siso}}; dB::Bool=false)
  Z, P = sys.num[1], sys.den[1]

  # To find w180, first solve imag(sys(jw)) == 0
  rZ, iZ, rP, iP = polysplit(Z)..., polysplit(P)...
  w180tmp = roots(iZ*rP - rZ*iP)

  # Evaluate the system at those frequencies
  gtmp = real(vec(freqresp(sys, real(w180tmp))[1]))

  # Find valid indexes: w180 must be real and positive, and the phase must be π
  i180 = (imag(w180tmp) .== 0) & (real(w180tmp) .>= 0) & (gtmp .< 0)
  !any(i180) && (return Inf, NaN)

  gm, w180 = abs(gtmp[i180]), real(w180tmp[i180])

  # Find the smallest ammount of gain that would make the system unstable
  igm = findmin(abs(log10(gm)))[2]
  return (dB ? -20*log10(gm[igm]) : 1/real(gm[igm])), w180[igm]
end
# For other representations, convert to transferfunction
gainmargin(sys::LtiSystem, dB::Bool=false) = gainmargin(tf(sys), dB=dB)


# Compute the real and imaginary parts of a polynomial assuming the argument is complex (=jw)
function polysplit{T<:Real}(p::Poly{T})
  rp = copy(coeffs(p))
  ip = copy(rp)

  # j^(even) is real, j^(odd) is complex
  rp[2:2:end] = 0
  ip[1:2:end] = 0

  # Make sure signs are right (j^3 = -j etc...)
  rp[3:4:end] *= -1
  ip[4:4:end] *= -1

  return Poly(rp), Poly(ip)
end
