subroutine space_group()
!=========================================================================
! find from the isogonal_point_group of lattice the space group symmetry operation
!=========================================================================
  use variables, only : atom_prim_crys, natom_prim, natom_prim_type, holo_ele, holo_ord, &
                        spa_ele, spa_ord, spa_axis, spa_symb, ntype, stdout
  use mat_trans, only : mat_det, mat_inv
  use constants, only : tpi, err_tol, err_tol3, err_tol33
  implicit none
!-------------------------------------------------------------------------
! local variables
  integer(kind=4) :: i, j, k, p, q, ia, ib, ic, ncc, itype, start
  real(kind=8) :: trace
  real(kind=8) :: atom_tem(3), vec(3), diff(3), spa_tem(3,4,holo_ord)
  character(len=2) :: symb_order(10)
  integer(kind=4),allocatable :: spa_axis_tem(:,:), cc_ord(:)
  real(kind=8),allocatable :: spa_ele_tem(:,:,:)
  character(len=2),allocatable :: spa_symb_tem(:), cc_symb(:)
  logical :: match(natom_prim)
  logical,allocatable :: is_used(:)
!-------------------------------------------------------------------------
  !
  do i = 1, holo_ord
    do j = 1, natom_prim_type(1)
      vec(:) = atom_prim_crys(:,j) - matmul(holo_ele(:,:,i),atom_prim_crys(:,1))
      ! make sure components in range 0~1
      vec(:) = vec(:) - dble(nint(vec(:)))
      do q = 1, 3
        if(vec(q) < -err_tol) vec(q) = vec(q) + 1.d0
      end do
      ! check if holo_ele(:,:,i) together with vec take each atom into its equivalence
      match = .false.
      do itype = 1, ntype
        if(itype == 1) then
          start = 1
        else
          start = sum(natom_prim_type(1:(itype-1)))+1
        end if
        do k = start, sum(natom_prim_type(1:itype))
          atom_tem(:) = matmul(holo_ele(:,:,i),atom_prim_crys(:,k)) + vec(:)
          do p = start, sum(natom_prim_type(1:itype))
            do ia = -3, 3
              do ib = -3, 3
                do ic = -3, 3
                  diff(:) = atom_tem(:) - atom_prim_crys(:,p) - dble((/ia,ib,ic/))
                  if(all(abs(diff) < err_tol3)) match(k) = .true.
                end do
              end do
            end do
          end do
        end do
      end do
      !
      ! if true, holo_ele(:,:,i) together with vec is one space symmetry operation
      if(all(match)) then
        spa_ord = spa_ord + 1
        spa_tem(1:3,1:3,spa_ord) = holo_ele(1:3,1:3,i)
        spa_tem(1:3,4,spa_ord) = vec(1:3)
        exit
      end if
      !
    end do ! j
  end do ! i
  !
  allocate(spa_ele(3,4,spa_ord))
  allocate(spa_axis(3,spa_ord))
  allocate(spa_symb(spa_ord))
  spa_ele(:,:,1:spa_ord) = spa_tem(:,:,1:spa_ord)
  !
  ! rotaion axis and symbol
  do i = 1, spa_ord
    trace = spa_ele(1,1,i) + spa_ele(2,2,i) + spa_ele(3,3,i)
    if(mat_det(spa_ele(:,1:3,i)) > 0.d0) then
      if(abs(trace-3.d0) < err_tol) then
        spa_symb(i) = 'E '
        spa_axis(:,i) = (/ 0, 0, 0 /)
      else if(abs(trace-2.d0) < err_tol) then
        spa_symb(i) = 'C6'
        call rot_axis(spa_ele(:,1:3,i),tpi/6.d0,spa_axis(:,i))
      else if(abs(trace-1.d0) < err_tol) then
        spa_symb(i) = 'C4'
        call rot_axis(spa_ele(:,1:3,i),tpi/4.d0,spa_axis(:,i))
      else if(abs(trace) < err_tol) then
        spa_symb(i) = 'C3'
        call rot_axis(spa_ele(:,1:3,i),tpi/3.d0,spa_axis(:,i))
      else if(abs(trace+1.d0) < err_tol) then
        spa_symb(i) = 'C2'
        call rot_axis(spa_ele(:,1:3,i),tpi/2.d0,spa_axis(:,i))
      end if
    else if(mat_det(spa_ele(:,1:3,i)) < 0.d0) then
      if(abs(trace+3.d0) < err_tol) then
        spa_symb(i) = 'I '
        spa_axis(:,i) = (/ 0, 0, 0 /)
      else if(abs(trace+2.d0) < err_tol) then
        spa_symb(i) = 'S3'
        call rot_axis(spa_ele(:,1:3,i),tpi/6.d0,spa_axis(:,i))
      else if(abs(trace+1.d0) < err_tol) then
        spa_symb(i) = 'S4'
        call rot_axis(spa_ele(:,1:3,i),tpi/4.d0,spa_axis(:,i))
      else if(abs(trace) < err_tol) then
        spa_symb(i) = 'S6'
        call rot_axis(spa_ele(:,1:3,i),tpi/3.d0,spa_axis(:,i))
      else if(abs(trace-1.d0) < err_tol) then
        spa_symb(i) = 'm '
        call rot_axis(spa_ele(:,1:3,i),tpi/2.d0,spa_axis(:,i))
      end if
    end if
  end do
  ! 
  ! rearrange the operations in order of E-C2-C3-C4-C6-I-m-S3-S4-S6
  symb_order = (/ 'E ','C2','C3','C4','C6','I ','m ','S6','S4','S3' /)
  allocate(spa_ele_tem(3,4,spa_ord), spa_symb_tem(spa_ord), spa_axis_tem(3,spa_ord))
  allocate(is_used(spa_ord), cc_symb(spa_ord), cc_ord(spa_ord))
  spa_ele_tem = spa_ele; spa_symb_tem = spa_symb; spa_axis_tem = spa_axis
  q = 0
  is_used(:) = .false.
  ncc = 0 ! number of conjugate class
  cc_ord(:) = 0 ! order of conjugate class
  do i = 1, 10
    do j = 1, spa_ord
      if(spa_symb_tem(j) == symb_order(i)) then
        if(is_used(j)) cycle
        ncc = ncc + 1
        cc_symb(ncc) = symb_order(i)
        do k = 1, spa_ord
          do p = 1, spa_ord
            if(all(abs(matmul(spa_ele_tem(:,1:3,k),matmul(spa_ele_tem(:,1:3,j),mat_inv(spa_ele_tem(:,1:3,k)))) &
            -spa_ele_tem(:,1:3,p)) < err_tol33) .and. (.NOT. is_used(p))) then
              is_used(p) = .true.
              cc_ord(ncc) = cc_ord(ncc) + 1
              q = q + 1
              spa_ele(:,:,q) = spa_ele_tem(:,:,p)
              spa_symb(q) = spa_symb_tem(p)
              spa_axis(:,q) = spa_axis_tem(:,p)
            end if
          end do
        end do
      end if
    end do
  end do
  !
  deallocate(spa_ele_tem,spa_symb_tem,spa_axis_tem)
  !
  ! rotation method to choose the axis direction of C2 and m
  ! C2
  do i = 2, spa_ord
    if(spa_symb(i) == 'C2') then
      do j = 1, i-1
        if(spa_symb(j) == 'C2') then
          do k = 1, 2
            if(all(dble(spa_axis(:,i))+dble(cshift(spa_axis(:,j),k)) < err_tol3)) spa_axis(:,i) = -spa_axis(:,i)
          end do
        end if
      end do
    end if
  end do
  !
  ! m
  if(spa_symb(spa_ord/2+1) == 'I ') then
    do i = spa_ord/2+1, spa_ord
      if(spa_symb(i) == 'm ') then
        do j = 1, spa_ord/2
          if(spa_symb(j) == 'C2' .and. all(spa_ele(:,1:3,i)+spa_ele(:,1:3,j) < err_tol33)) spa_axis(:,i) = spa_axis(:,j)
        end do
      end if
    end do
  else
    do i = 2, spa_ord
      if(spa_symb(i) == 'm ') then
        do j = 1, i-1
          if(spa_symb(j) == 'm ') then
            do k = 1, 2
              if(all(dble(spa_axis(:,i))+dble(cshift(spa_axis(:,j),k)) < err_tol3)) spa_axis(:,i) = -spa_axis(:,i)
            end do
          end if
        end do
      end if
    end do
  end if
  !
  ! print class-wise elements of space group
  write (stdout,'(a2,a,i2,a,i2,a)') cc_symb(1), '(',cc_ord(1),'-',cc_ord(1),'):'
  do j = 1, 3
    write (stdout,'(3i3,f8.3,a)') nint(spa_ele(1,1:3,1)), spa_ele(j,4,1), ' | '
  end do
  q = 1
  do i = 2, ncc
    q = q + cc_ord(i-1)
    write (stdout,'(a2,a,i2,a,i2,a)') cc_symb(i),'(',q,'-',q+cc_ord(i)-1,'):'
    do j = 1, 3
      do k = q, q+cc_ord(i)-1
        write (stdout,'(3i3,f8.3,a,$)') nint(spa_ele(j,1:3,k)), spa_ele(j,4,k), ' | '
      end do
      write (stdout,*)
    end do
  end do
  !
  deallocate(is_used, cc_symb, cc_ord)
  !
  return
end subroutine space_group
    
    
      
        
        
          