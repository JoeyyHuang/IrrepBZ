subroutine isogonal_pg(spa_ord,spa_ele,iso_pg)
!=========================================================================
! find the isogonal point group (point group of the space group) 
! based on the space group matrices
!=========================================================================
  use mat_trans, only : mat_det
  use constants, only : tol
  implicit none
!-------------------------------------------------------------------------
! arguments
  integer(kind=4),intent(in) :: spa_ord
  real(kind=8),intent(in) :: spa_ele(3,4,spa_ord)
  character(len=3),intent(out) :: iso_pg
! local variables
  integer(kind=4) :: i, j, k
  integer(kind=4) :: npro_list(7), nimp_list(7) ! number of proper/improper rotation in given space with trace -3,-2,-1,0,1,2,3
  real(kind=8) :: trace
  type tra
    character(len=3) :: symbol
    integer(kind=4) :: ntra_pro(7)   ! number of proper rotation with trace -3, -2, -1, 0, 1, 2, 3
    integer(kind=4) :: ntra_imp(7)   ! number of improper rotation with trace -3, -2, -1, 0, 1, 2, 3
  end type
  type(tra) :: point_group(32)
!-------------------------------------------------------------------------
  ! triclinic
  point_group(1)%symbol = 'C1'
  point_group(1)%ntra_pro = (/ 0, 0, 0, 0, 0, 0, 1 /); point_group(1)%ntra_imp = (/ 0, 0, 0, 0, 0, 0, 0 /)
  point_group(2)%symbol = 'Ci'
  point_group(2)%ntra_pro = (/ 0, 0, 0, 0, 0, 0, 1 /); point_group(2)%ntra_imp = (/ 1, 0, 0, 0, 0, 0, 0 /)
  ! monoclinic
  point_group(3)%symbol = 'C2'
  point_group(3)%ntra_pro = (/ 0, 0, 1, 0, 0, 0, 1 /); point_group(3)%ntra_imp = (/ 0, 0, 0, 0, 0, 0, 0 /)
  point_group(4)%symbol = 'Cs'
  point_group(4)%ntra_pro = (/ 0, 0, 0, 0, 0, 0, 1 /); point_group(4)%ntra_imp = (/ 0, 0, 0, 0, 1, 0, 0 /)
  point_group(5)%symbol = 'C2h'
  point_group(5)%ntra_pro = (/ 0, 0, 1, 0, 0, 0, 1 /); point_group(5)%ntra_imp = (/ 1, 0, 0, 0, 1, 0, 0 /)
  ! orthorhombic
  point_group(6)%symbol = 'C2v'
  point_group(6)%ntra_pro = (/ 0, 0, 1, 0, 0, 0, 1 /); point_group(6)%ntra_imp = (/ 0, 0, 0, 0, 2, 0, 0 /)
  point_group(7)%symbol = 'D2'
  point_group(7)%ntra_pro = (/ 0, 0, 3, 0, 0, 0, 1 /); point_group(7)%ntra_imp = (/ 0, 0, 0, 0, 0, 0, 0 /)
  point_group(8)%symbol = 'D2h'
  point_group(8)%ntra_pro = (/ 0, 0, 3, 0, 0, 0, 1 /); point_group(8)%ntra_imp = (/ 1, 0, 0, 0, 3, 0, 0 /)
  ! trigonal
  point_group(9)%symbol = 'C3'
  point_group(9)%ntra_pro = (/ 0, 0, 0, 2, 0, 0, 1 /); point_group(9)%ntra_imp = (/ 0, 0, 0, 0, 0, 0, 0 /)
  point_group(10)%symbol = 'S6'
  point_group(10)%ntra_pro = (/ 0, 0, 0, 2, 0, 0, 1 /); point_group(10)%ntra_imp = (/ 1, 0, 0, 2, 0, 0, 0 /)
  point_group(11)%symbol = 'C3v'
  point_group(11)%ntra_pro = (/ 0, 0, 0, 2, 0, 0, 1 /); point_group(11)%ntra_imp = (/ 0, 0, 0, 0, 3, 0, 0 /)
  point_group(12)%symbol = 'D3'
  point_group(12)%ntra_pro = (/ 0, 0, 3, 2, 0, 0, 1 /); point_group(12)%ntra_imp = (/ 0, 0, 0, 0, 0, 0, 0 /)
  point_group(13)%symbol = 'D3d'
  point_group(13)%ntra_pro = (/ 0, 0, 3, 2, 0, 0, 1 /); point_group(13)%ntra_imp = (/ 1, 0, 0, 2, 3, 0, 0 /)
  ! tetragonal
  point_group(14)%symbol = 'C4'
  point_group(14)%ntra_pro = (/ 0, 0, 1, 0, 2, 0, 1 /); point_group(14)%ntra_imp = (/ 0, 0, 0, 0, 0, 0, 0 /)
  point_group(15)%symbol = 'S4'
  point_group(15)%ntra_pro = (/ 0, 0, 1, 0, 0, 0, 1 /); point_group(15)%ntra_imp = (/ 0, 0, 2, 0, 0, 0, 0 /)
  point_group(16)%symbol = 'D2d'
  point_group(16)%ntra_pro = (/ 0, 0, 3, 0, 0, 0, 1 /); point_group(16)%ntra_imp = (/ 0, 0, 2, 0, 2, 0, 0 /)
  point_group(17)%symbol = 'C4h'
  point_group(17)%ntra_pro = (/ 0, 0, 1, 0, 2, 0, 1 /); point_group(17)%ntra_imp = (/ 1, 0, 2, 0, 1, 0, 0 /)
  point_group(18)%symbol = 'C4v'
  point_group(18)%ntra_pro = (/ 0, 0, 1, 0, 2, 0, 1 /); point_group(18)%ntra_imp = (/ 0, 0, 0, 0, 4, 0, 0 /)
  point_group(19)%symbol = 'D4'
  point_group(19)%ntra_pro = (/ 0, 0, 5, 0, 2, 0, 1 /); point_group(19)%ntra_imp = (/ 0, 0, 0, 0, 0, 0, 0 /)
  point_group(20)%symbol = 'D4h'
  point_group(20)%ntra_pro = (/ 0, 0, 5, 0, 2, 0, 1 /); point_group(20)%ntra_imp = (/ 1, 0, 2, 0, 5, 0, 0 /)
  ! hexagonal
  point_group(21)%symbol = 'C6'
  point_group(21)%ntra_pro = (/ 0, 0, 1, 2, 0, 2, 1 /); point_group(21)%ntra_imp = (/ 0, 0, 0, 0, 0, 0, 0 /)
  point_group(22)%symbol = 'C3h'
  point_group(22)%ntra_pro = (/ 0, 0, 0, 2, 0, 0, 1 /); point_group(22)%ntra_imp = (/ 0, 2, 0, 0, 1, 0, 0 /)
  point_group(23)%symbol = 'D3h'
  point_group(23)%ntra_pro = (/ 0, 0, 3, 2, 0, 0, 1 /); point_group(23)%ntra_imp = (/ 0, 2, 0, 0, 4, 0, 0 /)
  point_group(24)%symbol = 'C6h'
  point_group(24)%ntra_pro = (/ 0, 0, 1, 2, 0, 2, 1 /); point_group(24)%ntra_imp = (/ 1, 2, 0, 2, 1, 0, 0 /)
  point_group(25)%symbol = 'C6v'
  point_group(25)%ntra_pro = (/ 0, 0, 1, 2, 0, 2, 1 /); point_group(25)%ntra_imp = (/ 0, 0, 0, 0, 6, 0, 0 /)
  point_group(26)%symbol = 'D6'
  point_group(26)%ntra_pro = (/ 0, 0, 7, 2, 0, 2, 1 /); point_group(26)%ntra_imp = (/ 0, 0, 0, 0, 0, 0, 0 /)
  point_group(27)%symbol = 'D6h'
  point_group(27)%ntra_pro = (/ 0, 0, 7, 2, 0, 2, 1 /); point_group(27)%ntra_imp = (/ 1, 2, 0, 2, 7, 0, 0 /)
  ! cubic
  point_group(28)%symbol = 'T'
  point_group(28)%ntra_pro = (/ 0, 0, 3, 8, 0, 0, 1 /); point_group(28)%ntra_imp = (/ 0, 0, 0, 0, 0, 0, 0 /)
  point_group(29)%symbol = 'Th'
  point_group(29)%ntra_pro = (/ 0, 0, 3, 8, 0, 0, 1 /); point_group(29)%ntra_imp = (/ 1, 0, 0, 8, 3, 0, 0 /)
  point_group(30)%symbol = 'Td'
  point_group(30)%ntra_pro = (/ 0, 0, 3, 8, 0, 0, 1 /); point_group(30)%ntra_imp = (/ 0, 0, 6, 0, 6, 0, 0 /)
  point_group(31)%symbol = 'O'
  point_group(31)%ntra_pro = (/ 0, 0, 9, 8, 6, 0, 1 /); point_group(31)%ntra_imp = (/ 0, 0, 0, 0, 0, 0, 0 /)
  point_group(32)%symbol = 'Oh'
  point_group(32)%ntra_pro = (/ 0, 0, 9, 8, 6, 0, 1 /); point_group(32)%ntra_imp = (/ 1, 0, 6, 8, 9, 0, 0 /)
  !
  npro_list = 0
  nimp_list = 0
  do i = 1, spa_ord
    trace = spa_ele(1,1,i) + spa_ele(2,2,i) + spa_ele(3,3,i)
    if(mat_det(spa_ele(:,:,i)) > 0.d0) then
      do j = -3, 3
        if(abs(trace+dble(j)) < tol) npro_list(4-j) = npro_list(4-j) + 1
      end do
    else if(mat_det(spa_ele(:,:,i)) < 0.d0) then
      do k = -3, 3
        if(abs(trace+dble(k)) < tol) nimp_list(4-k) = nimp_list(4-k) + 1
      end do
    end if
  end do
  !
  do i = 1, 32
    if(all(npro_list == point_group(i)%ntra_pro) .and. all(nimp_list == point_group(i)%ntra_imp)) &
    iso_pg = point_group(i)%symbol
  end do
  !
  return
end subroutine isogonal_pg