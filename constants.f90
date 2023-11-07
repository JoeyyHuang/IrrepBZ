module constants
  implicit none
  public
  save
  real(kind=8),parameter :: tpi = 2.d0*acos(-1.d0)
  real(kind=8),parameter :: err_tol = 1.d-10   ! scalar error tolerance
  real(kind=8),parameter :: err_tol3(3) = 1.d-10   ! vector error tolerance
  real(kind=8),parameter :: err_tol22(2,2) = 1.d-10   ! 2*2 tensor error tolerance
  real(kind=8),parameter :: err_tol33(3,3) = 1.d-10   ! 3*3 tensor error tolerance
  complex(kind=8),parameter :: im = cmplx(0.d0,1.d0,8)
end module constants