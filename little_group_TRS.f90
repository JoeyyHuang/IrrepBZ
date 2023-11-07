subroutine little_order_t(k_crys, lit_ord, lit_ord_t)
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
  integer(kind=4),intent(out) :: lit_ord, lit_ord_t
! local variables
  integer(kind=4) :: i, ia, ib, ic
  real(kind=8) :: k_crys2(3)
!-------------------------------------------------------------------------
  ! number of symmetry operations that make k*alpha_i^-1 = \pm k + G_i
  lit_ord = 0; lit_ord_t = 0
  do i = 1, spa_ord
    k_crys2 = matmul(k_crys,mat_inv(spa_ele(:,1:3,i)))
    do ia = -3, 3
      do ib = -3, 3
        do ic = -3, 3
          if(all(abs(k_crys2-k_crys-(/ia,ib,ic/)) < err_tol3)) lit_ord = lit_ord + 1
          if(all(abs(k_crys2-k_crys-(/ia,ib,ic/)) < err_tol3)) lit_ord_t = lit_ord_t + 1
          if(all(abs(k_crys2+k_crys-(/ia,ib,ic/)) < err_tol3)) lit_ord_t = lit_ord_t + 1
        end do
      end do
    end do
  end do
  !
  return
end subroutine little_order_t
!
!
subroutine little_group_t(k_crys, lit_ord, lit_ord_t, lit_spa_t, lit_ele_t, lit_symb_t, lit_axis_t)
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
  integer(kind=4),intent(in) :: lit_ord, lit_ord_t
  integer(kind=4),intent(out) :: lit_axis_t(3,lit_ord_t), lit_spa_t(lit_ord_t)
  real(kind=8),intent(out) :: lit_ele_t(3,4,lit_ord_t)
  character(len=2),intent(out) :: lit_symb_t(lit_ord_t)
! local variables
  integer(kind=4) :: i, m, n, ia, ib, ic
  real(kind=8) :: k_crys2(3)
!-------------------------------------------------------------------------
  m = 0; n = 0
  ! symmetry operations that make k*alpha_i^-1 = k + G_i
  do i = 1, spa_ord
    k_crys2 = matmul(k_crys,mat_inv(spa_ele(:,1:3,i)))
    do ia = -3, 3
      do ib = -3, 3
        do ic = -3, 3
          if(all(abs(k_crys2-k_crys-(/ia,ib,ic/)) < err_tol3)) then
            m = m + 1
            lit_spa_t(m) = i
            lit_ele_t(:,:,m) = spa_ele(:,:,i)
            lit_symb_t(m) = spa_symb(i)
            lit_axis_t(:,m) = spa_axis(:,i)
          end if
        end do
      end do
    end do
  end do
  n = m
  !
  ! symmetry operations that make k*alpha_i^-1 = -k + G_i
  do i = 1, spa_ord
    k_crys2 = matmul(k_crys,mat_inv(spa_ele(:,1:3,i)))
    do ia = -3, 3
      do ib = -3, 3
        do ic = -3, 3
          if(all(abs(k_crys2+k_crys-(/ia,ib,ic/)) < err_tol3)) then
            n = n + 1
            lit_ele_t(:,:,n) = spa_ele(:,:,i)
            lit_symb_t(n) = spa_symb(i)
            lit_axis_t(:,n) = spa_axis(:,i)
          end if
        end do
      end do
    end do
  end do
  !
  if(m /= lit_ord .or. n /= lit_ord_t) then
    write (stdout,*) 'Not consistent order of little cogroup'
    stop
  end if
  !
  return
end subroutine little_group_t