"""
`care(A, B, Q, R)`

Compute 'X', the solution to the continuous-time algebraic Riccati equation,
defined as A'X + XA - (XB)R^-1(B'X) + Q = 0, where R is non-singular.

Algorithm taken from:
Laub, "A Schur Method for Solving Algebraic Riccati Equations."
http://dspace.mit.edu/bitstream/handle/1721.1/1301/R-0859-05666488.pdf
"""
# TO DO: Change code by SLICOT version
function care{T<:BlasFloat}(A::StridedMatrix{T}, B::StridedMatrix{T},
            Q::StridedMatrix{T}, R::StridedMatrix{T})
    G = try
        B*inv(R)*B'
    catch
        error("R must be non-singular.")
    end

    Z = [A  -G;
        -Q  -A']

    S = schurfact(Z)
    S = ordschur(S, S.values.<0)
    U = S.Z

    (m, n) = size(U)
    U11 = U[1:div(m, 2), 1:div(n,2)]
    U21 = U[div(m,2)+1:m, 1:div(n,2)]
    return U21/U11
end
care{T1<:Real,T2<:Real,T3<:Real,T4<:Real}(A::StridedMatrix{T1}, B::StridedMatrix{T2},
  Q::StridedMatrix{T3}, R::StridedMatrix{T4}) = care(float(A), float(B), float(Q), float(R))

"""
`dare(A, B, Q, R)`

Compute `X`, the solution to the discrete-time algebraic Riccati equation,
defined as A'XA - X - (A'XB)(B'XB + R)^-1(B'XA) + Q = 0, where A and R
are non-singular.

Algorithm taken from:
Laub, "A Schur Method for Solving Algebraic Riccati Equations."
http://dspace.mit.edu/bitstream/handle/1721.1/1301/R-0859-05666488.pdf
"""
# TO DO: Change code by SLICOT version
function dare{T<:BlasFloat}(A::StridedMatrix{T}, B::StridedMatrix{T},
            Q::StridedMatrix{T}, R::StridedMatrix{T})
    G = try
        B*inv(R)*B'
    catch
        error("R must be non-singular.")
    end

    Ait = try
        inv(A)'
    catch
        error("A must be non-singular.")
    end

    Z = [A + G*Ait*Q   -G*Ait;
         -Ait*Q        Ait]

    S = schurfact(Z)
    S = ordschur(S, S.values.*S.values.<=1)
    U = S.Z

    (m, n) = size(U)
    U11 = U[1:div(m, 2), 1:div(n,2)]
    U21 = U[div(m,2)+1:m, 1:div(n,2)]
    return U21/U11
end
dare{T1<:Real,T2<:Real,T3<:Real,T4<:Real}(A::StridedMatrix{T1}, B::StridedMatrix{T2},
  Q::StridedMatrix{T3}, R::StridedMatrix{T4}) = dare(float(A), float(B), float(Q), float(R))
