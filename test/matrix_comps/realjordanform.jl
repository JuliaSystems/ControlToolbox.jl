# test real jordan form
D = [-1, -1+2im, -1+1im, -1-1im,-1-2im, -1+2im, -1-2im]
C1 = [-1 -1; 1 -1]
C2 = [-1 -2; 2 -1]
JD = zeros(Int, 7,7)
JD[1,1] = -1
JD[2:3,2:3] = C1
JD[4:5,4:5] = JD[6:7,6:7] = C2
JD[4:5,6:7] = eye(Int, 2)

@test realjordanform(D) == JD
