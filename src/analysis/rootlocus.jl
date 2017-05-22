type RootLocusResponse{T} <: SystemResponse
  K::Vector{T}        # gains
  systf::TransferFunction{Val{:siso}}
  real_p::Matrix{T}   # real part of poles
  imag_p::Matrix{T}   # imaginary part of poles (one row for each k in K)
  real_m0::T          # center
  imag_m0::T          # center
  d0::T               # radius
  tloc::Int

  function (::Type{RootLocusResponse}){S<:Real,U<:Real,V<:Real}(
    K::Vector{S}, systf::TransferFunction{Val{:siso}}, real_p::Matrix{U},
    imag_p::Matrix{V}, real_m0::Real, imag_m0::Real, d0::Real)
    if size(real_p) != size(imag_p)
      warn("RootLocusResponse: mag and phase must have same dimensions")
      throw(DomainError())
    end
    if length(K) != size(real_p,1)
      warn("RootLocusResponse: size(real_p,1) ≠ length(K)")
      throw(DomainError())
    end

    T = promote_type(S,U,V,typeof(real_m0),typeof(imag_m0), typeof(d0))
    new{T}(Vector{T}(K),
           systf,
           Matrix{T}(real_p),
           Matrix{T}(imag_p),
           convert(T,real_m0),
           convert(T,imag_m0),
           convert(T,d0),
           0)
  end

  @compat function (rls::RootLocusResponse)(k::Real)
    r = rls.systf.mat[1]
    nump = num(r)
    denp = den(r)
    roots(denp + k*nump)
  end
end

# Iteration interface
start(rls::RootLocusResponse)       = (rls.tloc = 1)
done(rls::RootLocusResponse, state) = state >= length(rls.K)
next(rls::RootLocusResponse, state) = (state+=1; rls.tloc = state; (rls, state))

# Plot some outputs for some inputs
@recipe function f(rls::RootLocusResponse)

  r = rls.systf.mat[1]
  lastidx = length(rls.K)
  endidx  = rls.tloc == 0 ? lastidx : rls.tloc
  startp  = roots(den(r))
  endp    = roots(num(r))

  # Define plotting rules
  layout                :=  (1,1)
  legend                := :none

  xlabel                --> "Real Axis"
  ylabel                --> "Imaginary Axis"
  xlims                 --> (min(minimum(rls.real_p), rls.real_m0-1.1rls.d0),
                             max(maximum(rls.real_p), rls.real_m0+1.1rls.d0))
  ylims                 --> (min(minimum(rls.imag_p), rls.imag_m0-1.1rls.d0),
                             max(maximum(rls.imag_p), rls.imag_m0+1.1rls.d0))
  title                 --> "Root locus"
  legend                --> :topright
  grid                  --> true

  # series for root locus line
  for idx in indices(rls.real_p, 2)
    @series begin
      subplot           :=  1
      hover             --> [@sprintf "(K = %.2f, p = %.2f%+.2f im)" rls.K[k] rls.real_p[k,idx] rls.imag_p[k,idx] for k in eachindex(rls.K)]

      rls.real_p[1:endidx,idx], rls.imag_p[1:endidx,idx]
    end
  end

  # start pole locations
  @series begin
    subplot             :=  1
    markershape         := :xcross
    linealpha           := 0
    markerstrokecolor   := :black
    markercolor         := :black
    markerstrokewidth   := 2
    markersize          := 6

    real(startp), imag(startp)
  end

  # endpoint circles at lastidx
  @series begin
    subplot             :=  1
    markershape         := :circle
    linealpha           := 0
    markerstrokecolor   := :black
    markercolor         := :white
    markersize          := 6

    real(endp), imag(endp)
  end
end

"""
    `rootlocus(sys, [K])` -> rl

Computes the root locus response `rl` of `sys`.

Returns an object that contains the poles of the closed-loop system
`feedback(sys, k)`, for `k` in `K`. The function `rootlocus` can also be
called with

    `rootlocus(nump::Poly, denp::Poly, [K])` -> rl

`nump` and `denp` are the numerator and denominator polynomials
respectively for the open-loop system.

`lr` is a custom data type (`RootLocusResponse`) containing

  * Gains (lr.K) over which the root locus is taken
  * System (lr.systf) such that the root locus corresponds to the poles of
    feedback(lr.systf, k)
  * Real part of poles (lr.real_p)
  * Imaginary part of poles (lr.real_p)

Plotting recipe is defined for `lr`, which allows for
`plot(lr; <keyword arguments>)` when `using Plots`.

# Examples
```julia
julia> sys = tf([2., 5, 1], [1., 2, 3]);

julia> rl = rootlocus(sys);
"""
function rootlocus{S}(sys::LtiSystem{Val{:siso}}, K::AbstractVector{S}=Float64[];
  kwargs...)
  systf = tf(sys)
  rootlocus(systf, K; kwargs...)
end

function rootlocus{S<:Real}(systf::TransferFunction{Val{:siso}},
  K::AbstractVector{S}=Float64[]; N::Int=100)
  r  = systf.mat[1]
  nump = num(r)
  denp = den(r)
  Nr = 100
  m0 = mean(vcat(roots(nump), roots(denp)))
  d0 = max(maxabs(zeros(r)-m0), maxabs(poles(r)-m0))
  if isempty(K)
    k = 1e-4
    while maxabs(roots(denp + k*nump)) < 10d0
      k *= 10
      if k > 1e2
        break
      end
    end
    N = max(N, 100Nr)
    K = logspace(-3, log10(k), N)
  else
    N = length(K)
    Nr = min(Nr, N)
  end

  C = complex(promote_type(eltype(r)..., Float16, eltype(K)))
  # find the rootlocus for N points
  nbroots   = maximum(degree(r))
  prevroots = roots(denp)
  pvec      = zeros(C, N, nbroots)
  for i in eachindex(K)
    k = K[i]
    currroots = roots(denp + k*nump)
    nextroots = copy(currroots)
    for j in eachindex(prevroots)
      p        = prevroots[j]
      _,l      = findmin(norm.(currroots-p))
      pk       = splice!(currroots,l)
      pvec[i,j] = pk
    end
    prevroots = nextroots
  end

  # find index set for N good points of the rootlocus
  k_inds = ones(Int,Nr)
  totdiff = sum(norm(pvec[i,:]-pvec[i-1,:]) for i in 2:N)
  currdiff = 0.0
  k = 1
  for i in 2:N
    currdiff += norm(pvec[i,:]-pvec[i-1,:])
    if (currdiff > k/Nr*totdiff) && (k < Nr-1)
      k += 1
      k_inds[k] = i
    end
  end
  k_inds[Nr] = N

  return RootLocusResponse(K[k_inds],
                           systf,
                           real(pvec[k_inds,:]),
                           imag(pvec[k_inds,:]),
                           real(m0),
                           imag(m0),
                           d0)
end
