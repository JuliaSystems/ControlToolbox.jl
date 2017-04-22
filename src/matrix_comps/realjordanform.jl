"""
    realjordanform(D) -> J

Returns the real Jordan canonical form of the eigenvalues specified in
the vector D.

# Examples
```julia
julia> D = [-1, -1+2im, -1-2im, -1+2im, -1-2im];
julia> realjordanform(D)
5x5 Array{Int64,2}:
 -1   0   0   0   0
  0  -1  -2   1   0
  0   2  -1   0   1
  0   0   0  -1  -2
  0   0   0   2  -1
```
"""

function realjordanform{T<:Number}(D::AbstractVector{T})
  sort!(D, by=imag)
  sort!(D, by=x->abs(imag(x))) # make sure conjugated pairs are next to eachother
  sort!(D, by=real)

  J = zeros(real(T), length(D), length(D))
  i = 1
  while i <= length(D)
    d = D[i]
    if isreal(d)
      n = findlast(D, d) - i + 1 # multiplicity of real eigenvalue
      @assert all(D[i+j] == d for j = 1:n-1) "There should be $n eigenvalues equal to $d"

      # fill jordan block for real eigenvalue
      J[i,i] = d
      for k = 1:n-1
        J[i+k-1,i+k] = 1
        J[i+k,i+k]   = d
      end
      i += n
    else
      n = findlast(D, d) - i + 1 # multiplicity of complex eigenvalue
      @assert all(D[i+j] == d for j = 1:n-1) "There should be $n eigenvalues equal to $d"
      @assert all(D[i+n+j] == conj(d) for j = 1:n-1) "There should be $n eigenvalues equal to $(conj(d))"
      # construct real block C
      a = real(d)
      b = imag(d)
      C = [a b; -b a]

      # fill jordan block for complex eigenvalue
      J[i:i+1,i:i+1] = C
      for k = 1:n-1
        J[i+2(k-1)+(0:1),i+2k+(0:1)] = eye(C)
        J[i+2k+(0:1),i+2k+(0:1)]     = C
      end
      i += 2n
    end
  end
  J
end
