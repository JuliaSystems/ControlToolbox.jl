"""
    observability(sys) -> O

Compute the observability matrix `O` for `sys`, the `StateSpace` system description.
"""
function observability(sys::StateSpace)
  nx      = numstates(sys)
  ny, nu  = size(sys)
  temp    = zeros(nx*ny, nx)
  temp[1:ny, :] = sys.C
  for i = 1:nx-1
    temp[(1 + i*ny):(1 + i)*ny, :] = temp[((i - 1)*ny + 1):i*ny, :] * sys.A
  end
  temp
end
