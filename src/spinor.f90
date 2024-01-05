subroutine spinor(lit_ord, lit_symb, lit_axis, lit_ele, fac_spin)
!=========================================================================
! obtain the factor system stemming from the spinor representation of little group
!=========================================================================
  use variables, only : a_prim
  use constants, only : im, tpi, tol, tol22, tol33, stdout
  implicit none
!-------------------------------------------------------------------------
! arguments
  integer(kind=4),intent(in) :: lit_ord
  integer(kind=4),intent(in) :: lit_axis(3,lit_ord)
  real(kind=8),intent(in) :: lit_ele(3,4,lit_ord)
  character(len=2),intent(in) :: lit_symb(lit_ord)
  real(kind=8),intent(out) :: fac_spin(lit_ord,lit_ord)
! local variables
  integer(kind=4) :: i, j, k
  real(kind=8) :: angle
  complex(kind=8) :: SU2_mat(2,2,lit_ord), vec(3)
!-------------------------------------------------------------------------
! SU2 associated with SO3, in Pauli gauge SU2(R) = SU2(IR)
  SU2_mat = (0.d0,0.d0)
  do i = 1, lit_ord
    if(lit_symb(i) == 'E ' .or. lit_symb(i) == 'I ') then
      SU2_mat(:,:,i) = reshape((/(1.d0,0.d0), (0.d0,0.d0), (0.d0,0.d0), (1.d0,0.d0)/),(/2, 2/))
    else
      if(lit_symb(i) == 'C2' .or. lit_symb(i) == 'm ') angle = tpi/2.d0
      if(lit_symb(i) == 'C3' .or. lit_symb(i) == 'S6') angle = tpi/3.d0
      if(lit_symb(i) == 'C4' .or. lit_symb(i) == 'S4') angle = tpi/4.d0
      if(lit_symb(i) == 'C6' .or. lit_symb(i) == 'S3') angle = tpi/6.d0
      vec = matmul(a_prim,dble(lit_axis(:,i)))
      vec = vec / sqrt(vec(1)**2+vec(2)**2+vec(3)**2)
      SU2_mat(:,:,i) = reshape((/cos(angle/2.d0)-im*vec(3)*sin(angle/2.d0), (vec(2)-im*vec(1))*sin(angle/2.d0), &
                      -(vec(2)+im*vec(1))*sin(angle/2.d0), cos(angle/2.d0)+im*vec(3)*sin(angle/2.d0)/),(/2, 2/))
    end if
  end do
  ! factor system in spinor projective representation
  fac_spin = 0.d0
  do i = 1, lit_ord
    do j = 1, lit_ord
      do k = 1, lit_ord
        if(all(abs(matmul(lit_ele(:,1:3,i),lit_ele(:,1:3,j))-lit_ele(:,1:3,k)) < tol33)) then
          if(all(abs(matmul(SU2_mat(:,:,i),SU2_mat(:,:,j))-SU2_mat(:,:,k)) < tol22)) fac_spin(i,j) = 1.d0
          if(all(abs(matmul(SU2_mat(:,:,i),SU2_mat(:,:,j))+SU2_mat(:,:,k)) < tol22)) fac_spin(i,j) = -1.d0
          if(abs(abs(fac_spin(i,j))-1.d0) > tol) then
            write (stdout,*) 'Something wrong in SU2 representation'
            stop
          end if
        end if
      end do
    end do
  end do
  !
  return
end subroutine spinor