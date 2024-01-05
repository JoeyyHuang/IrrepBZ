module read_binary
!=========================================================================
! read the wavefunction in binary format obtained from QE calculation.
! Adapted from Modules/io_base.f90 of Quantum Espresso
!=========================================================================
  contains
  subroutine read_a_wfc(filename, dims, evc)
  implicit none 
  !
  character (len=*), intent(in)  :: filename 
  integer(kind=4), intent(out), allocatable :: dims(:,:)   
  complex(kind=8), intent(out), allocatable :: evc(:,:)
  ! 
  integer :: iuni = 1111, i, nbnd, ispin, npol, ngw, igwx, ik
  real(kind=8) :: scalef
  real(kind=8) :: b1(3), b2(3), b3(3), xk(3)
  logical :: gamma_only 
  !
  open(UNIT = iuni, FILE = trim(filename), FORM = 'unformatted', status = 'old') 
  read(iuni) ik, xk, ispin, gamma_only, scalef
  read(iuni) ngw, igwx, npol, nbnd
  read(iuni) b1, b2, b3 
  !
  allocate (dims(3,igwx))
  read (iuni) dims (1:3,1:igwx)
  !
  allocate (evc(npol*igwx,nbnd))
  do i = 1, nbnd 
    read(iuni) evc (1:npol*igwx,i) 
  end do 
  close(iuni) 
  end subroutine read_a_wfc 
end module read_binary