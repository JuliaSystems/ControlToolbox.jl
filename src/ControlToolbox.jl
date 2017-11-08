module ControlToolbox

using Compat
using Optim
using Polynomials
using RationalFunctions
using RecipesBase
using LTISystems

import Base: step, norm
import Base.LinAlg: BlasFloat
import Base: start, next, done
import LTISystems: LtiSystem, StateSpace, TransferFunction, SystemResponse
import MathProgBase: eval_grad_f, eval_f, eval_g, features_available, initialize

using MathProgBase: loadproblem!, NonlinearModel, setwarmstart!, optimize!
using MathProgBase: status, AbstractNLPEvaluator, AbstractMathProgSolver

export
  c2d,
  Discretization,
  isstable,
  rootlocus,
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
  gram,
  realjordanform,
  Poleplacement,
  place,
  controllability,
  observability,
  dcgain

include("c2d.jl")
include("d2c.jl")
include("analysis/controllability.jl")
include("analysis/damp.jl")
include("analysis/dcgain.jl")
include("analysis/isstable.jl")
include("analysis/margins.jl")
include("analysis/markovparam.jl")
include("analysis/observability.jl")
include("analysis/rootlocus.jl")
include("matrix_comps/riccati.jl")
include("matrix_comps/lyapunov.jl")
include("matrix_comps/covar.jl")
include("matrix_comps/norm.jl")
include("matrix_comps/gram.jl")
include("matrix_comps/realjordanform.jl")
include("design/place.jl")

end # module
