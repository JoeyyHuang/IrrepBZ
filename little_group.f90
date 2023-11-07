subroutine little_order(k_crys, lit_ord)
!=========================================================================
! get the order of little group of a given k-point within FBZ
!=========================================================================
  use variables, only : spa_ord, spa_ele
  use mat_trans, only : mat_inv
  use constants, only : err_tol3
  implicit none
!-------------------------------------------------------------------------
! arguments
  real(kind=8),intent(in) :: k_crys(3)
  integer(kind=4),intent(out) :: lit_ord
! local variables
  integer(kind=4) :: i, ia, ib, ic
  real(kind=8) :: k_crys2(3)
!-------------------------------------------------------------------------
  ! number of symmetry operations that make k*alpha_i^-1 = k + G_i
  lit_ord = 0
  do i = 1, spa_ord
    k_crys2 = matmul(k_crys,mat_inv(spa_ele(:,1:3,i)))
    do ia = -3, 3
      do ib = -3, 3
        do ic = -3, 3
          if(all(abs(k_crys2-k_crys-(/ia,ib,ic/)) < err_tol3)) lit_ord = lit_ord + 1
        end do
      end do
    end do
  end do
  !
  return
end subroutine little_order
!
!
subroutine little_group(k_crys, lit_ord, lit_spa, lit_ele, lit_symb, lit_axis, lit_cg)
!=========================================================================
! get the little co-group of a given k-point within FBZ
!=========================================================================
  use variables, only : stdout, spa_ord, spa_ele, spa_symb, spa_axis
  use mat_trans, only : mat_inv
  use constants, only : err_tol3
  implicit none
!-------------------------------------------------------------------------
! arguments
  real(kind=8),intent(in) :: k_crys(3)
  integer(kind=4),intent(in) :: lit_ord
  integer(kind=4),intent(out) :: lit_axis(3,lit_ord), lit_spa(lit_ord)
  real(kind=8),intent(out) :: lit_ele(3,4,lit_ord)
  character(len=2),intent(out) :: lit_symb(lit_ord)
  character(len=3),intent(out) :: lit_cg
! local variables
  integer(kind=4) :: i, n, ia, ib, ic
  real(kind=8) :: k_crys2(3)
!-------------------------------------------------------------------------
  ! symmetry operations that make k*alpha_i^-1 = k + G_i
  n = 0
  do i = 1, spa_ord
    k_crys2 = matmul(k_crys,mat_inv(spa_ele(:,1:3,i)))
    do ia = -3, 3
      do ib = -3, 3
        do ic = -3, 3
          if(all(abs(k_crys2-k_crys-(/ia,ib,ic/)) < err_tol3)) then
            n = n + 1
            lit_spa(n) = i
            lit_ele(:,:,n) = spa_ele(:,:,i)
            lit_symb(n) = spa_symb(i)
            lit_axis(:,n) = spa_axis(:,i)
          end if
        end do
      end do
    end do
  end do
  !
  if(n /= lit_ord) then
    write (stdout,*) 'Not consistent order of little cogroup'
    stop
  end if
  !
  call isogonal_pg(lit_ord,lit_ele,lit_cg)
  !
  return
end subroutine little_group