module ControlToolbox

using ControlCore
using Polynomials

import Base: step
import Base.LinAlg: BlasFloat
import ControlCore: LtiSystem, StateSpace, RationalTF, ZeroPoleGain
import ControlCore: Siso, Continuous

export
  c2d,
  lsim,
  step,
  impulse,
  isstable,
  rlocus,
  damp,
  dampreport,
  phasemargin,
  gainmargin,
  markovparam,
  care,
  dare,
  clyap,
  dlyap,
  covar,
  norm,
  gram

# using DSP

include("c2d.jl")
include("d2c.jl")
include("analysis/isstable.jl")
include("analysis/margins.jl")
#include("analysis/rlocus.jl")
include("analysis/damp.jl")
include("analysis/dcgain.jl")
include("analysis/markovparam.jl")
include("matrix_comps/riccati.jl")
include("matrix_comps/lyapunov.jl")
include("matrix_comps/covar.jl")
include("matrix_comps/norm.jl")
include("matrix_comps/gram.jl")
#include("simulation/utils.jl")
#include("simulation/lsim.jl")
#include("simulation/step.jl")
#include("simulation/impulse.jl")

end # module
