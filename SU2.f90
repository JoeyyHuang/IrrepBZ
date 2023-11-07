subroutine SU2(spa_ord, spa_symb, spa_axis, spa_ele, fac_spin)
!=========================================================================
! SU2 representation of the symmetry group in Pauli gauge
!=========================================================================
  use variables, only : stdout, a_prim
  use constants, only : im, tpi, err_tol, err_tol22, err_tol33
  implicit none
!-------------------------------------------------------------------------
! arguments
  integer(kind=4),intent(in) :: spa_ord
  integer(kind=4),intent(in) :: spa_axis(3,spa_ord)
  real(kind=8),intent(in) :: spa_ele(3,4,spa_ord)
  character(len=2),intent(in) :: spa_symb(spa_ord)
  real(kind=8),intent(out) :: fac_spin(spa_ord,spa_ord)
! local variables
  integer(kind=4) :: i, j, k
  real(kind=8) :: angle
  complex(kind=8) :: SU2_mat(2,2,spa_ord), vec(3)
!-------------------------------------------------------------------------
! SU2 associated with SO3, in Pauli gauge SU2(R) = SU2(IR)
  SU2_mat = (0.d0,0.d0)
  do i = 1, spa_ord
    if(spa_symb(i) == 'E ' .or. spa_symb(i) == 'I ') then
      SU2_mat(:,:,i) = reshape((/(1.d0,0.d0), (0.d0,0.d0), (0.d0,0.d0), (1.d0,0.d0)/),(/2, 2/))
    else
      if(spa_symb(i) == 'C2' .or. spa_symb(i) == 'm ') angle = tpi/2.d0
      if(spa_symb(i) == 'C3' .or. spa_symb(i) == 'S6') angle = tpi/3.d0
      if(spa_symb(i) == 'C4' .or. spa_symb(i) == 'S4') angle = tpi/4.d0
      if(spa_symb(i) == 'C6' .or. spa_symb(i) == 'S3') angle = tpi/6.d0
      vec = matmul(a_prim,dble(spa_axis(:,i)))
      vec = vec / sqrt(vec(1)**2+vec(2)**2+vec(3)**2)
      SU2_mat(:,:,i) = reshape((/cos(angle/2.d0)-im*vec(3)*sin(angle/2.d0), (vec(2)-im*vec(1))*sin(angle/2.d0), &
                      -(vec(2)+im*vec(1))*sin(angle/2.d0), cos(angle/2.d0)+im*vec(3)*sin(angle/2.d0)/),(/2, 2/))
    end if
  end do
  ! factor system in spinor projective representation
  fac_spin = 0.d0
  do i = 1, spa_ord
    do j = 1, spa_ord
      do k = 1, spa_ord
        if(all(abs(matmul(spa_ele(:,1:3,i),spa_ele(:,1:3,j))-spa_ele(:,1:3,k)) < err_tol33)) then
          if(all(abs(matmul(SU2_mat(:,:,i),SU2_mat(:,:,j))-SU2_mat(:,:,k)) < err_tol22)) fac_spin(i,j) = 1.d0
          if(all(abs(matmul(SU2_mat(:,:,i),SU2_mat(:,:,j))+SU2_mat(:,:,k)) < err_tol22)) fac_spin(i,j) = -1.d0
          if(abs(abs(fac_spin(i,j))-1.d0) > err_tol) then
            write (stdout,*) 'Something wrong in SU2 representation'
            stop
          end if
        end if
      end do
    end do
  end do
  !
  return
end subroutine SU2