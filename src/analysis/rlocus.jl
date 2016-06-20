"""
`rlocus(sys, [K])`
Computes the root locus of `sys`.

Returns the closed-loop system
- poles at `K`, if `K` is a number
- pole trajectories, if `K` is a vector of gains
- pole trajectories at evenly spaced intervals of the gain, if `K` is
omitted
"""
function rlocus(sys::ControlCore.LtiSystem)
  N::Integer = 100    # Number of points in root locus

  K = Vector{Real}(N) # Vector of gains
  fill!(K, 1e8)
  K[1] = 1e-4         # Initial gain

  α = 1.05            # Maximum increase of pole distance
  λ = 2.0             # Decrease in gain per iteration

  last_max = max_cplx_distance(rlocus(sys, K[1])[1])
  for n in 2:N
      max_d = max_cplx_distance(rlocus(sys, K[n])[1])
      while max_d >= α * last_max
          K[n] /= λ
          max_d = max_cplx_distance(rlocus(sys, K[n])[1])
      end
      last_max = max_d
  end

  return rlocus(sys, K)
end

function rlocus{T<:Real}(sys::ControlCore.LtiSystem, K::AbstractVector{T})
  plist = Matrix{eltype(poles(sys))}(length(K), length(denvec(sys))-1)
  for idx in eachindex(K)
    plist[idx,:] = rlocus(sys, K[idx])[1]
  end
  (plist, K)
end

rlocus(sys::ControlCore.LtiSystem, k::Real) = poles(feedback(sys,k)), k

"""
`max_cplx_distance(C)`
Computes the maximal distance between any two complex points in the set `C`.
"""
function max_cplx_distance{T<:Complex}(C::AbstractVector{T})
    max_d = 0.0
    for p in C
      for q in C
          if abs(p-q) >= max_d
              max_d = abs(p-q)
          end
      end
    end
    return max_d
end
