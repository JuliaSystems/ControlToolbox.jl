include("../src/ControlToolbox.jl")

A = [[.5,.1] [.1,.5]]
B = reshape([2,1],2,1)
C = reshape([1,0],1,2)
D = zeros(Float64,1,1)

s = ss(A,B,C,D,1)
u = randn(10,1)

s2 = ss(A,B,C,D)
s2d = c2d(s2, 0.1)[1]
s2d = c2d(s2, 0.1, Discretization.BackwardEuler())[1]
s2d = c2d(s2, 0.1, Discretization.ForwardEuler())[1]
s2d = c2d(s2, 0.1, Discretization.GeneralizedBilinear(0.3))[1]
s2d = c2d(s2, 0.1, Discretization.FOH())[1]

t = Array(1:length(u))
(y,x) = lsim(s,u,t)
(y,x) = lsim(s2,u,t)

s2 = tf([1,2,1],[2,1,1],1)
x0 = zeros(Float64, max(length(denvec(s2)), length(numvec(s2)))-1)
y2 = lsim(s2, randn(10,1), t, x0)

s3 = zpk([1,2,3],[2,5,7],6.7,1)

u2 = randn(20,1)
t2 = Array(1:length(u2))
y3 = lsim(s3, u2, t2)



s4 = mimo(reshape([s3,s3,s2,s2, s2,s2,s3,s3],2,4))
lsim(s4, randn(20,4), t2)

s = ss(A,B,C,D,1)
s2 = ss(A,B,C,D)
step(s,t)
step(s2,t)

impulse(s,t)
impulse(s2,t)


s2 = tf([1,2,1],[2,1,1])
s3 = zpk([1,2,3],[2,5,7],6.7)
s4 = mimo(reshape([s3,s3,s2,s2, s2,s2,s3,s3],2,4))
lsim(s4, randn(20,4), t2)


s3 = zpk([1,2,3],[0.1,-0.2,0.4],6.7,1)
isstable(s3)
impulse(s3,t)
y, t, x = impulse(s)

t0 = 0:0.05:2
systf = mimo([tf([1],[1,1]) 0; 0 tf([1,-1],[1,2,1])])
sysss = ss([-1 0 0; 0 -2 -1; 0 1 0], [1 0; 0 1; 0 0], [1 0 0; 0 1 -1], zeros(2,2))
yreal = zeros(length(t0), 2, 2)
xreal = zeros(length(t0), 3, 2)

#Step tf
ds = c2d(systf, 0.1)

y, t, x = step(systf, Array(t0))
yreal[:,1,1] = 1-exp(-t)
yreal[:,2,2] = -1+exp(-t)+2*exp(-t).*t
typeof(systf) <: ControlCore.LtiSystem

nu = numinputs(systf)
ny = numoutputs(systf)

y, t, x = step(sysss, Array(t0))

norm(y[:,1,1]-yreal[:,1,1])

y[:,1,1]

yreal[:,1,1]
