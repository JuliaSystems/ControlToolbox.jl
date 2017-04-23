module ControlToolbox

using SystemsBase
using Polynomials
using Optim

import Base: step
import Base.LinAlg: BlasFloat
import SystemsBase: LtiSystem, StateSpace, RationalTF

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
  gram,
  realjordanform,
  place

include("c2d.jl")
include("d2c.jl")
include("analysis/isstable.jl")
include("analysis/margins.jl")
include("analysis/rlocus.jl")
include("analysis/damp.jl")
include("analysis/dcgain.jl")
include("analysis/markovparam.jl")
include("matrix_comps/riccati.jl")
include("matrix_comps/lyapunov.jl")
include("matrix_comps/covar.jl")
include("matrix_comps/norm.jl")
include("matrix_comps/gram.jl")
include("matrix_comps/realjordanform.jl")
include("simulation/utils.jl")
#include("simulation/lsim.jl")
#include("simulation/step.jl")
#include("simulation/impulse.jl")
include("design/place.jl")

end # module
