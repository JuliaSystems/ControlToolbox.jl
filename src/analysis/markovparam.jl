"""
`markovparam(sys, n)`

Compute the `n`th markov parameter of state-space system `sys`. This is defined
as the following:

`h(0) = D`

`h(n) = C*A^(n-1)*B`
"""
function markovparam(sys::StateSpace, n::Integer)
    n < 0 && error("n must be >= 0")
    return n == 0 ? sys.D : sys.C * sys.A^(n-1) * sys.B
end

markovparam(sys::LtiSystem, n::Integer) = markovparam(ss(sys),n)
