"""
    place(A, B, p[; rtol = 1e-3, α::Real = 0, iter::Integer=200]) -> K

Places the closed-loop poles of the single- or multi-input system

x = Ax + Bu

at the poles of the vector `p` of desired self-conjugate closed-loop pole
locations. Computes a gain matrix `K` such that the state feedback `u = –Kx`
places the closed-loop poles at the locations `p`. In other words, the
eigenvalues of `A + BK` match the entries of `p` (up to the ordering).
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
function place{S<:Real,U<:Real,V<:Number}(
  A::AbstractArray{S, 2}, B::AbstractArray{U, 2}, p::AbstractArray{V, 1};
  α::Real = 0.99, rtol::Real=1e-3, iter::Integer=200)

  # define optimization problem
  d = Poleplacement(A, B, p; α=α, rtol=rtol)

  if(size(A,1) - d.m  > 0) # if there is something to optimize over
    df = OnceDifferentiable(x -> eval_f(d::Poleplacement, x),
                           (x, grad_f) -> eval_grad_f(d::Poleplacement, grad_f, x))
    res = optimize(df,
                  rand(eltype(d.xₚ), length(d.xₚ)),
                  Optim.BFGS(),
                  Optim.Options(iterations = iter))
  end
  # transform controller according to the factorization of A
  K = hcat(zeros(size(d.K,1), d.m), d.K)*d.Qₐ.'
end

function place{S<:Real,U<:Real,V<:Number}(
  A::AbstractArray{S, 2}, B::AbstractArray{U, 2}, p::AbstractArray{V, 1},
  solver::AbstractMathProgSolver;
  α::Real = 0.99, rtol::Real=1e-3, iter::Integer=200)

  # define optimization problem
  d = Poleplacement(A, B, p; α=α, rtol=rtol)

  if(size(A,1) - d.m  > 0) # if there is something to optimize over
    G0 = randn(size(A,1)-d.m, size(B,2)) # random initial point
    n = length(G0)
    l = -Inf*ones(n)
    u = Inf*ones(n)
    lb = Float64[]
    ub = Float64[]
    numconst = 0

    m = NonlinearModel(solver)
    loadproblem!(m, n, numconst, l, u, lb, ub, :Min, d)
    setwarmstart!(m, vec(G0))
    optimize!(m)
    if status(m) != :Optimal
      warn("place: solution not optimal")
      throw(InvalidStateException())
    end
  end
  # transform controller according to the factorization of A
  K = hcat(zeros(size(d.K,1), d.m), d.K)*d.Qₐ.'
end

type Poleplacement{T,M} <: AbstractNLPEvaluator
  K::M
  A::M
  B::M
  Ã::M
  α::T
  xₚ::Vector{T}
  X::M
  Xᵢ::M
  m::Integer
  Qₐ::M

  @compat function (::Type{Poleplacement}){S<:Real,U<:Real,V<:Number}(
    A::AbstractArray{S, 2}, B::AbstractArray{U, 2}, p::AbstractArray{V, 1};
    α::Real = 0.99, rtol::Real=1e-3)
    T = promote_type(eltype(A), eltype(B), real(eltype(p)), typeof(α), Float16)
    A = one(T)*A
    B = one(T)*B
    p = one(T)*p
    α = one(T)*α

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
    p = p[!select]
    Ã = realjordanform(p)

    # move eigenvalues to be moved to the lower right corner of the real schur
    # form, i.e. in A₃
    # QₐᵀAQₐ = Rₐ =  [A₁ A₂;      QₐᵀB = QB = [B₁;
    #                 0  A₃]                   B₂]
    ordschur!(Rₐ, Qₐ, select)
    Qᵦ = Qₐ.'*B
    A₃ = Rₐ[m+1:end,m+1:end]
    B₂ = Qᵦ[m+1:end,:]

    # minimize Frobenius norm of state feedback K, ‖K‖ and condition number
    # of X, κ(X) according to cost function J = ακ(X) + (1-α)‖K‖²
    G  = randn(T,size(B₂,2),size(A₃,2))
    xₚ = randn(T,length(G))
    X  = eye(T,size(A₃,1))
    Xᵢ = eye(T,size(A₃,1))
    K  = zeros(T,G*Xᵢ)

    new{T,typeof(X)}(K, A₃, B₂, Ã, α, xₚ, X, Xᵢ, m, Qₐ)
  end
end

# helper function to save work done in both function call and gradient call
function common_first_sylvester!(d::Poleplacement, x)
  if x != d.xₚ
    d.xₚ[:] = x
    G = reshape(x, size(d.B,2), size(d.A,2))
    # solve AX - XÃ = -B*G
    X, scale = LAPACK.trsyl!('N', 'N', d.A, d.Ã, -d.B*G, -1)
    scale!(X, inv(scale)) # this scale is always 1?
    d.X[:]  = X
    d.Xᵢ[:] = inv(X)
    d.K[:]  = G*d.Xᵢ
  end
end

function initialize(d::Poleplacement, requested_features::Vector{Symbol})
    for feat in requested_features
        if !(feat in [:Grad])
            error("Unsupported feature $feat")
        end
    end
end

features_available(d::Poleplacement) = [:Grad]

function eval_f(d::Poleplacement, x)
  common_first_sylvester!(d, x)
  d.α*(sumabs2(d.X) + sumabs2(d.Xᵢ))/2 + (1-d.α)*sumabs2(d.K)/2
end

eval_g(d::Poleplacement, g, x) = nothing

function eval_grad_f(d::Poleplacement, grad_f, x)
  common_first_sylvester!(d, x)
  X,Xᵢ,K,α = d.X,d.Xᵢ,d.K,d.α
  # solve ÃU - UA = -S
  H = Xᵢ*K'
  S = α*(-X.' + Xᵢ*Xᵢ.'*Xᵢ) + (1-α)*H*K  # Xi*Xi.'*Xi (X.'*X)\Xi
  U,scale = LAPACK.trsyl!('N', 'N', d.Ã, d.A, -S, -1)
  scale!(U, inv(scale)) # this scale is always 1?

  # calculate gradient
  grad_f[:] = vec((1-α)*H' + d.B.'*U.')
end
