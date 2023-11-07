module mat_trans
  implicit none
  !
  contains
function mat_inv(mat)
!=========================================================================
! get the inversion matrix based on direct-reciprocal lattice relationship
!=========================================================================
  implicit none
!-------------------------------------------------------------------------
! arguments
  real(kind=8),intent(in) :: mat(3,3)
  real(kind=8) :: mat_inv(3,3)
! local variables
  real(kind=8) :: omega
!-------------------------------------------------------------------------
  ! volume (a1xa2).a3
  omega = (mat(1,2)*mat(2,3)-mat(2,2)*mat(1,3))*mat(3,1) - &
          (mat(1,1)*mat(2,3)-mat(2,1)*mat(1,3))*mat(3,2) + &
          (mat(1,1)*mat(2,2)-mat(2,1)*mat(1,2))*mat(3,3)
  ! inversion matrix
  mat_inv(1,1) = (mat(2,2)*mat(3,3)-mat(3,2)*mat(2,3))/omega
  mat_inv(2,1) = (mat(3,1)*mat(2,3)-mat(2,1)*mat(3,3))/omega
  mat_inv(3,1) = (mat(2,1)*mat(3,2)-mat(3,1)*mat(2,2))/omega
  mat_inv(1,2) = (mat(3,2)*mat(1,3)-mat(1,2)*mat(3,3))/omega
  mat_inv(2,2) = (mat(1,1)*mat(3,3)-mat(3,1)*mat(1,3))/omega
  mat_inv(3,2) = (mat(3,1)*mat(1,2)-mat(1,1)*mat(3,2))/omega
  mat_inv(1,3) = (mat(1,2)*mat(2,3)-mat(2,2)*mat(1,3))/omega
  mat_inv(2,3) = (mat(2,1)*mat(1,3)-mat(1,1)*mat(2,3))/omega
  mat_inv(3,3) = (mat(1,1)*mat(2,2)-mat(2,1)*mat(1,2))/omega
  !
  return
end function mat_inv
!
function mat_det(mat)
!=========================================================================
! get the matrix determinant
!=========================================================================
  implicit none
!-------------------------------------------------------------------------
! arguments
  real(kind=8),intent(in) :: mat(3,3)
  real(kind=8) :: mat_det
!-------------------------------------------------------------------------
  mat_det = (mat(2,2)*mat(3,3)-mat(3,2)*mat(2,3))*mat(1,1) - &
            (mat(2,1)*mat(3,3)-mat(3,1)*mat(2,3))*mat(1,2) + &
            (mat(2,1)*mat(3,2)-mat(3,1)*mat(2,2))*mat(1,3)
 
  return
end function mat_det
!
end module mat_trans