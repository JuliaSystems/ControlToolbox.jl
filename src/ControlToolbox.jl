module ControlToolbox

using ControlCore

import Base: step
import ControlCore: LtiSystem, StateSpace, Siso, Continuous

export c2d, lsim, step, impulse, isstable, rlocus

# using DSP

include("c2d.jl")
#include("simulation/utils.jl")
#include("simulation/lsim.jl")
#include("simulation/step.jl")
#include("simulation/impulse.jl")
#include("analysis/isstable.jl")
#include("analysis/rlocus.jl")
#include("analysis/dcgain.jl")

end # module
