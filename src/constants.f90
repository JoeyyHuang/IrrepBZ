module constants
  implicit none
  public
  save
  real(kind=8),parameter :: pi = acos(-1.d0)
  real(kind=8),parameter :: tpi = 2.d0*acos(-1.d0)
  real(kind=8),parameter :: tol = 1.d-5   ! scalar error tolerance
  real(kind=8),parameter :: tol3(3) = 1.d-5   ! vector error tolerance
  real(kind=8),parameter :: tol22(2,2) = 1.d-5   ! 2*2 tensor error tolerance
  real(kind=8),parameter :: tol33(3,3) = 1.d-5   ! 3*3 tensor error tolerance
  complex(kind=8),parameter :: im = cmplx(0.d0,1.d0,8)
  real(kind=8),parameter :: bohr2ang = 0.52917721067d0
  real(kind=8),parameter :: c = 0.262465831d0 ! 2m/hbar**2 in units of 1/eV Ang^2
  ! output file units
  integer(kind=4),parameter :: stdout = 100
  integer(kind=4),parameter :: qe_band = 200
  integer(kind=4),parameter :: qe_phonon = 300
  integer(kind=4),parameter :: vasp_band = 400
  integer(kind=4),parameter :: vasp_phonon = 500
end module constants