subroutine little_order(k_crys, lit_ord)
!=========================================================================
! get the order of little co-group at a given k-point within FBZ
!=========================================================================
  use variables, only : spa_ord, spa_ele
  use mat_trans, only : mat_inv
  use constants, only : tol3
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
          if(all(abs(k_crys2-k_crys-(/ia,ib,ic/)) < tol3)) lit_ord = lit_ord + 1
        end do
      end do
    end do
  end do
  !
  return
end subroutine little_order
!
!
subroutine little_group(funit, k_crys, lit_ord, lit_spa, lit_ele, lit_symb, lit_axis)
!=========================================================================
! get the little co-group at a given k-point within FBZ
!=========================================================================
  use variables, only : spa_ord, spa_ele, spa_symb, spa_axis
  use mat_trans, only : mat_inv
  use constants, only : tol3, tol33
  implicit none
!-------------------------------------------------------------------------
! arguments
  real(kind=8),intent(in) :: k_crys(3)
  integer(kind=4),intent(in) :: funit, lit_ord
  integer(kind=4),intent(out) :: lit_axis(3,lit_ord), lit_spa(lit_ord)
  real(kind=8),intent(out) :: lit_ele(3,4,lit_ord)
  character(len=2),intent(out) :: lit_symb(lit_ord)
! local variables
  integer(kind=4) :: i, j, k, p, q, ncc, ia, ib, ic
  integer(kind=4) :: lit_axis_cc(3,lit_ord), lit_spa_cc(lit_ord), cc_ord(lit_ord)
  real(kind=8) :: k_crys2(3), lit_ele_cc(3,4,lit_ord)
  character(len=2) :: symb_order(10), lit_symb_cc(lit_ord), cc_symb(lit_ord)
  character(len=3) :: lit_cg
  logical:: is_used(lit_ord)
!-------------------------------------------------------------------------
  ! symmetry operations that make k*alpha_i^-1 = k + G_i
  q = 0
  do i = 1, spa_ord
    k_crys2 = matmul(k_crys,mat_inv(spa_ele(:,1:3,i)))
    do ia = -3, 3
      do ib = -3, 3
        do ic = -3, 3
          if(all(abs(k_crys2-k_crys-(/ia,ib,ic/)) < tol3)) then
            q = q + 1
            lit_spa(q) = i
            lit_ele(:,:,q) = spa_ele(:,:,i)
            lit_symb(q) = spa_symb(i)
            lit_axis(:,q) = spa_axis(:,i)
          end if
        end do
      end do
    end do
  end do
  !
  if(q /= lit_ord) then
    write (funit,*) 'Error: not consistent order of little co-group'
    stop
  end if
  !
  ! conjugate class
  symb_order = (/ 'E ','C2','C3','C4','C6','I ','m ','S6','S4','S3' /)
  q = 0
  is_used(:) = .false.
  ncc = 0 ! number of conjugate class
  cc_ord(:) = 0 ! order of conjugate class
  do i = 1, 10
    do j = 1, lit_ord
      if(lit_symb(j) == symb_order(i)) then
        if(is_used(j)) cycle
        ncc = ncc + 1
        cc_symb(ncc) = symb_order(i)
        do p = 1, lit_ord
          do k = 1, lit_ord
            if(all(abs(matmul(lit_ele(:,1:3,k),matmul(lit_ele(:,1:3,j),mat_inv(lit_ele(:,1:3,k)))) &
            -lit_ele(:,1:3,p)) < tol33) .and. (.NOT. is_used(p))) then
              is_used(p) = .true.
              cc_ord(ncc) = cc_ord(ncc) + 1
              q = q + 1
              lit_spa_cc(q) = lit_spa(p)
              lit_ele_cc(:,:,q) = lit_ele(:,:,p)
              lit_symb_cc(q) = lit_symb(p)
              lit_axis_cc(:,q) = lit_axis(:,p)
              exit
            end if
          end do
        end do
      end if
    end do
  end do
  !
  ! the little co-group
  call isogonal_pg(lit_ord,lit_ele,lit_cg)
  if(funit == 300 .or. funit == 500) then
    write (funit,'(a,a)') 'The little co-group at q-point: ', lit_cg
  else
    write (funit,'(a,a)') 'The little co-group at k-point: ', lit_cg
  end if
  !
  ! print class-wise elements of little group
  write (funit,'(a,i3)') 'Order of little co-group:', lit_ord
  write (funit,'(a,i3)') 'Number of conjugate class:', ncc
  write (funit,*)
  write (funit,'(a)') 'The representative operations of little group in each conjugate class:'
  write (funit,'(a)') 'numbers in square bracket indicate the labels of elements in space group'
  write (funit,'(a2,a,i2,a)') cc_symb(1),'[',lit_spa_cc(1),']:'
  do j = 1, 3
    write (funit,'(3i3,2x,i3,f8.3,a)') nint(lit_ele_cc(j,1:3,1)), lit_axis(j,1), lit_ele_cc(j,4,1)+1.d-8, ' |'
  end do
  q = 1
  do i = 2, ncc
    q = q + cc_ord(i-1)
    write (funit,'(a2,a,$)') cc_symb(i),'['
    if(q == q+cc_ord(i)-1) then
      write (funit,'(i2,$)') lit_spa_cc(q)
    else
      do k = q, q+cc_ord(i)-2
        write (funit,'(i2,a,$)') lit_spa_cc(k),','
      end do
      write (funit,'(i2,$)') lit_spa_cc(q+cc_ord(i)-1)
    end if
    write (funit,'(a)') ']:'
    do j = 1, 3
      do k = q, q+cc_ord(i)-1
        write (funit,'(3i3,2x,i3,f8.3,a,$)') nint(lit_ele_cc(j,1:3,k)), lit_axis(j,k), lit_ele_cc(j,4,k)+1.d-8, ' |'
      end do
      write (funit,*)
    end do
  end do
  write (funit,*)
  !
  return
end subroutine little_group