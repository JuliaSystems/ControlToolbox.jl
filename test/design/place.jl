# These examples are mostly taken from [1] and the references therein.
# References
#
#  [1]: Chu, E.K. "Pole assignment via the Schur form." Systems & control
#       letters 56.4 (2007): 303-314.
# example 1
A₁ = [1.38       -2.077e-1  6.715  -5.67   ;
     -5.814e-1   -4.29      0       6.77e-1;
      1.067       4.273    -6.654   5.893  ;
      4.8e-2      4.273     1.343  -2.104  ]
B₁ = [0      0    ;
      5.679  0    ;
      1.136 -3.146;
      1.136  0    ]
p₁ = [-0.2, -0.5, -5.0566, -8.6659]
κ₁ = 5
c₁ = 5

# example 2
A₂ = [-1.094e-1   6.28e-2    0         0         0       ;
       1.306     -2.132      9.87e-1   0         0       ;
       0          1.595     -3.149     1.547     0       ;
       0          3.55e-2    2.632    -4.257     1.855   ;
       0          2.3e-3     0         1.636e-1 -1.625e-1]
B₂ = [0          0       ;
      6.38e-2    0       ;
      8.38e-2   -1.396e-1;
      1.004e-1  -2.06e-1 ;
      6.3e-3    -1.28e-2 ]
p₂ = [-0.2, -0.5, -1, -1+im, -1-im]
κ₂ = 1_000
c₂ = 200

# example 3
A₃ = [-65   65  -19.5  19.4;
       0.1 -0.1  0     0   ;
       1    0   -0.5  -1   ;
       0    0    0.4   0   ]
B₃ = [65  0  ;
      0   0  ;
      0   0.4;
      0   0  ]
p₃ = [-1, -2, -3, -4.0]
κ₃ = 50
c₃ = 100

# example 4
A₄ = [ 0   1   0;
       0   0   1;
      -6  -11 -6]
B₄ = [1 1;
      0 1;
      1 1]
p₄ = [-1, -2, -3]
κ₄ = 1_000
c₄ = 1

# example 5
D  = diagm([1,10,0.1,0.1,10])
D₁ = diagm([1,0.1,10,10,0.1])
D₂ = diagm([1e5,1e2])
A₅ = D*[-1.29e-1  0  3.96e-2   2.5e-2    1.91e-2;
         3.29e-3  0 -7.79e-5   1.22e-4  -6.21e-1;
         7.18e-2  0 -1e-1      8.87e-4  -3.85   ;
         4.11e-2  0  0        -8.22e-2   0      ;
         3.51e-4  0  3.5e-5    4.26e-5  -7.43e-2]*D₁
B₅ = D*[0        1.39e-3;
        0        3.59e-5;
        0       -9.89e-3;
        2.49e-5  0      ;
        0       -5.34e-6]*D₂
p₅ = [-0.01, -0.02, -0.03, -0.04, -0.05]
κ₅ = 1_000
c₅ = 5

# example 6
A₆ = [5.8765 9.3456 4.5634 9.3520;
      6.6526 0.5867 3.5829 0.6534;
      0      9.6738 7.4876 4.7654;
      0      0      6.6784 2.5678]
B₆ = [3.9878 0.5432;
      0      2.7650;
      0      0     ;
      0      0     ]
p₆ = [-29.4986, -10.0922, 2.5201+6.891im, 2.5201-6.891im]
κ₆ = 5
c₆ = 50

# example 7 (pole of multiplicity 4 > rank(B) = 2)
A₇ = [1.38       -2.077e-1  6.715  -5.67   ;
     -5.814e-1   -4.29      0       6.77e-1;
      1.067       4.273    -6.654   5.893  ;
      4.8e-2      4.273     1.343  -2.104  ]
B₇ = [0      0    ;
      5.679  0    ;
      1.136 -3.146;
      1.136  0    ]
p₇ = [-0.5+1im, -0.5-1im, -0.5+1im, -0.5-1im]
κ₇ = 1e12
c₇ = 5

# example 8
A₈ = [ 0        1        0       0      ;
       5.32e-7 -4.18e-1 -0.12    2.32e-2;
      -4.62e-9  1       -0.752  -2.39e-2;
      -5.61e-1  0       -0.3    -1.74e-2]
B₈ = [ 0         0      ;
      -1.72e-1   7.45e-6;
      -2.38e-2  -7.78e-5;
       0         3.69e-3]
p₈ = [-1,-2,-3,-4]
κ₈ = 1_000
c₈ = 1_000

# example 9 (ill-conditioned example)
A₉ = [0 0  0   0   ;
      1 10 100 1000;
      0 1  10  100 ;
      0 0  1   10  ]
B₉ = [1 0;
      0 1;
      0 0;
      0 0]
p₉ = [-1,-2,-3,-4]
κ₉ = 1e6
c₉ = 10_000

for (A,B,p,κ,c) in ((A₁,B₁,p₁,κ₁,c₁), (A₂,B₂,p₂,κ₂,c₂), (A₃,B₃,p₃,κ₃,c₃),
                    (A₄,B₄,p₄,κ₄,c₄), (A₅,B₅,p₅,κ₅,c₅), (A₆,B₆,p₆,κ₆,c₆),
                    (A₇,B₇,p₇,κ₇,c₇), (A₈,B₈,p₈,κ₈,c₈), (A₉,B₉,p₉,κ₉,c₉))
  F = place(A, B, p, solver=NLoptSolver(algorithm=:LD_MMA))
  F = place(A, B, p)
  λ, V = eig(A+B*F)

  @test cond(V)     < κ
  @test vecnorm(F)  < c
  for i in 1:length(p)
    j = findfirst(y->isapprox(y, p[i], rtol=1e-3), λ)
    if j > 0
      splice!(λ,j)
    end
  end
  @test isempty(λ)
end
