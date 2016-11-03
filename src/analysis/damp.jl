"""
`Wn, zeta, ps = damp(sys)`

Compute the natural frequencies, `Wn`, and damping ratios, `zeta`, of the
poles, `ps`, of `sys`
"""
function damp(sys::LtiSystem)
  ps = poles(sys)
  if !iscontinuous(sys)
    #Ts = sys.Ts == -1 ? 1 : sys.Ts
    ps = log(ps)/sys.Ts
  end
  Wn = abs(ps)
  order = sortperm(Wn)
  Wn = Wn[order]
  ps = ps[order]
  zeta = -cos(angle(ps))
  return Wn, zeta, ps
end

"""
`dampreport(sys)`

Display a report of the poles, damping ratio, natural frequency, and time
constant of the system `sys`
"""
# TO DO: @printf currently throws an error if there are complex poles...
function dampreport(io::IO, sys::LtiSystem)
    Wn, zeta, ps = damp(sys)
    t_const = 1./(Wn.*zeta)
    header =
    ("|     Pole      |   Damping     |   Frequency   | Time Constant |\n"*
     "|               |    Ratio      |   (rad/sec)   |     (sec)     |\n"*
     "+---------------+---------------+---------------+---------------+")
    println(io, header)
    for i=1:length(ps)
        p, z, w, t = ps[i], zeta[i], Wn[i], t_const[i]
        @printf(io, "|  %-13.3e|  %-13.3e|  %-13.3e|  %-13.3e|\n", p, z, w, t)
    end
end

dampreport(sys::LtiSystem) = dampreport(STDOUT, sys)
