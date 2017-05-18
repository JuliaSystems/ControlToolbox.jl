"""
    controllability(sys) -> C

Compute the controllability matrix `C` for `sys`, the `StateSpace` system description.
"""
function controllability(sys::StateSpace)
  nx      = numstates(sys)
  ny, nu  = size(sys)
  temp    = zeros(nx, nx*nu)
  temp[:, 1:nu] = sys.B
  for i = 1:nx-1
    temp[:, (1 + i*nu):(1 + i)*nu] = sys.A * temp[:, ((i - 1)*nu + 1):i*nu]
  end
  temp
end
