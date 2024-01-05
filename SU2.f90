subroutine SU2()
!=========================================================================
! SU2 representation of the symmetry group in Pauli gauge
!=========================================================================
  use variables, only : spa_ord, spa_axis, spa_ele, spa_symb, a_prim
  use mat_trans, only : mat_inv
  use constants, only : im, tpi, tol, tol22, tol33, stdout
  implicit none
!-------------------------------------------------------------------------
! local variables
  integer(kind=4) :: i, j, k, p, q, ncc
  integer(kind=4),allocatable :: cc_ord(:)
  real(kind=8) :: angle
  complex(kind=8) :: SU2_mat(2,2,spa_ord), vec(3)
  character(len=2) :: symb_order(10)
  character(len=2),allocatable :: cc_symb(:)
  logical,allocatable :: is_used(:)
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
  !
  symb_order = (/ 'E ','C2','C3','C4','C6','I ','m ','S6','S4','S3' /)
  allocate(is_used(spa_ord), cc_symb(spa_ord), cc_ord(spa_ord))
  is_used(:) = .false.
  ncc = 0 ! number of conjugate class
  cc_ord(:) = 0 ! order of conjugate class
  do i = 1, 10
    do j = 1, spa_ord
      if(spa_symb(j) == symb_order(i)) then
        if(is_used(j)) cycle
        ncc = ncc + 1
        cc_symb(ncc) = symb_order(i)
        do k = 1, spa_ord
          do p = 1, spa_ord
            if(all(abs(matmul(spa_ele(:,1:3,k),matmul(spa_ele(:,1:3,j),mat_inv(spa_ele(:,1:3,k)))) &
            -spa_ele(:,1:3,p)) < tol33) .and. (.NOT. is_used(p))) then
              is_used(p) = .true.
              cc_ord(ncc) = cc_ord(ncc) + 1
            end if
          end do
        end do
      end if
    end do
  end do
  ! print class-wise elements of space group and their SU2 matrxi representations
  write (stdout,'(a)') 'The SU2 matrix representations in Pauli gauge: SU2(R) = SU2(IR)'
  write (stdout,'(a)') 'The factor system due to spinor in the projective representation is determined by these matrices'
  write (stdout,'(a2,a,i2,a,i2,a)') cc_symb(1), '(',cc_ord(1),'-',cc_ord(1),'):'
  do j = 1, 2
    write (stdout,'(2("(",f6.3,","f6.3,") "),a)') SU2_mat(j,:,1)+(1.d-8,1.d-8), '| '
  end do
  q = 1
  do i = 2, ncc
    q = q + cc_ord(i-1)
    write (stdout,'(a2,a,i2,a,i2,a)') cc_symb(i),'(',q,'-',q+cc_ord(i)-1,'):'
    do j = 1, 2
      do k = q, q+cc_ord(i)-1
        write (stdout,'(2("(",f6.3,","f6.3,") "),a,$)') SU2_mat(j,:,k)+(1.d-8,1.d-8), '| '
      end do
      write (stdout,*)
    end do
  end do
  write (stdout,*)
  deallocate(is_used, cc_symb, cc_ord)
  !
  return
end subroutine SU2