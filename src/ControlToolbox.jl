module ControlToolbox

using ControlCore
using Polynomials

import Base: step
import ControlCore: LtiSystem, StateSpace, RationalTF, Siso, Continuous

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
  gainmargin

# using DSP

include("c2d.jl")
include("d2c.jl")
#include("simulation/utils.jl")
#include("simulation/lsim.jl")
#include("simulation/step.jl")
#include("simulation/impulse.jl")
include("analysis/isstable.jl")
include("analysis/margins.jl")
#include("analysis/rlocus.jl")
include("analysis/damp.jl")
include("analysis/dcgain.jl")
include("analysis/markovparam.jl")

end # module
