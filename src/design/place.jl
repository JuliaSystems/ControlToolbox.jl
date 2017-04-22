"""
    place(A, B, p[; rtol = 1e-3, α::Real = 0, iter::Integer=200]) -> K

Places the closed-loop poles of the single- or multi-input system

x = Ax + Bu

at the poles of the vector `p` of desired self-conjugate closed-loop pole
locations. Computes a gain matrix `K` such that the state feedback `u = –Kx`
places the closed-loop poles at the locations `p`. In other words, the
eigenvalues of `A – BK` match the entries of `p` (up to the ordering).
Place works for multi-input systems and is based on the algorithm from [1].
A closed loop pole `pₖ` in `p` that is close enough to an eigenvalue of `A`,
`λₖ` will not be moved by the feedback. The tolerance is set `|λₖ-pₖ|< rtol`.

The non-uniqueness of solutions is used to minimize the sensitivity of
closed-loop eigenvalues and the norm of the corresponding state feedback matrix
`K`. The parameter `α ∈ [0,1]` determines the relative weight of these
criteria. `α=0` only minimizes the norm of `K` and `α=1` only minimizes the
sensitivity of the closed loop eigenvalues measured by the condition number of
the eigenvectors of the closed loop eigenvalues.

You can also use place for estimator gain selection by transposing the `A` matrix
and substituting `C.'` for `B`.

l = place(A.',C.',p).'

# Examples
```julia
julia> A = [0 0  0   0   ;
            1 10 100 1000;
            0 1  10  100 ;
            0 0  1   10  ];

julia> B = [1 0;
            0 1;
            0 0;
            0 0];

julia> p = [-1, -2, -1+1im, -1-1im];

julia> K = place(A, B, p);

julia> λ, V = eig(A+B*K);

julia> λ
Complex{Float64}[4]
  -2.00 + 0.00im
  -1.00… + 1.00…im
  -1.00… - 1.00…im
  -1.00… + 0.00im
```

# References

-  [1]: Varga, A. "Robust pole assignment via Sylvester equation based state
        feedback parametrization." IEEE International Symposium on
        Computer-Aided Control System Design, 2000.
"""
function place{M1<:AbstractMatrix,M2<:AbstractMatrix,M3<:AbstractVector}(
    A::M1, B::M2, p::M3; α::Real = 0.99, rtol::Real=1e-3, iter::Integer=200)

  # /TODO check (A,B) controllable
  if length(p) != size(A,1)
    warn("The number of poles in `p` does not match match the size of `A`")
    throw(DomainError())
  end

  # schur factorization of A
  Rₐ, Qₐ, λ = schur(A)

  # find eigenvalues not moved
  select = zeros(Bool, length(p))
  for i in eachindex(p)
    j = findfirst(y->isapprox(p[i],y; rtol=rtol), λ)
    if (j > 0)
      splice!(λ, j)
      select[i] = true
    end
  end

  m = sum(select) # number of eigenvalues not moved
  p = float(p[!select]) # convert to float (A and B are handled by schur)
  Ã = realjordanform(p)

  # move eigenvalues to be moved to the lower right corner of the real schur
  # form, i.e. in A₃
  # QₐᵀAQₐ = Rₐ =  [A₁ A₂;      QₐᵀB = QB = [B₁;
  #                 0  A₃]                   B₂]
  ordschur!(Rₐ, Qₐ, select)
  QB = Qₐ.'*B
  A₃ = Rₐ[m+1:end,m+1:end]
  B₂ = QB[m+1:end,:]

  # find good G that parametrizes the freedom of the state feedback
  G = randn(eltype(A₃), size(B₂,2), size(A₃,2))
  # /TODO check (Ã,G) is observable

  # minimize Frobenius norm of state feedback K, ‖K‖ and condition number
  # of X, κ(X) according to cost function J = ακ(X) + (1-α)‖K‖²
  last_g = zeros(length(G))
  X      = zeros(A₃)
  Xi     = zeros(A₃)
  K      = zeros(G*Xi)

  df = OnceDifferentiable(x -> place_cost(x, last_g, X, Xi, K, A₃, Ã, B₂, α),
                         (x, stor) -> place_gradient!(stor, x, last_g, X, Xi, K, A₃, Ã, B₂, α))
  res = optimize(df,
                vec(G),
                Optim.BFGS(),
                Optim.Options(iterations = iter))

  # transform controller according to the factorization of A
  K = hcat(zeros(size(K,1), m), K)*Qₐ.'
  K
end

function place_cost(x, last_x, X, Xi, K, A, Ã, B, α)
  common_first_sylvester!(x, last_x, X, Xi, K, A, Ã, B)
  α*(sumabs2(X) + sumabs2(Xi))/2 + (1-α)*sumabs2(K)/2
end

function place_gradient!(storage, x, last_x, X, Xi, K, A, Ã, B, α)
  common_first_sylvester!(x, last_x, X, Xi, K, A, Ã, B)

  # solve ÃU - UA = -S
  H = Xi*K'
  S = α*(-X.' + Xi*Xi.'*Xi) + (1-α)*H*K  # Xi*Xi.'*Xi (X.'*X)\Xi
  U,scale = LAPACK.trsyl!('N', 'N', Ã, A, -S, -1)
  scale!(U, inv(scale)) # this scale is always 1?

  # calculate gradient
  storage[1:length(x)] = vec((1-α)*H' + B.'*U.')
end

# helper function to save work done in both function call and gradient call
function common_first_sylvester!(x, last_x, X, Xi, K, A, Ã, B)
  if x != last_x
    copy!(last_x, x)
    G = reshape(x, size(B,2), size(A,2))

    # solve AX - XÃ = -B*G
    X[:], scale = LAPACK.trsyl!('N', 'N', A, Ã, -B*G, -1)
    scale!(X, inv(scale)) # this scale is always 1?
    Xi[:] = inv(X)
    K[:]  = G*Xi
  end
end
