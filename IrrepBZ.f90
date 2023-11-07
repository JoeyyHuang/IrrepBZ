program IrrepBZ
!=========================================================================
! input: poscar file
! outut: the irreducible presentations and character table of little group at
! any given k-points in BZ, the induced irreducible presentations and character
! table of space group, both with and without considering the spin and time 
! reversal symmetry.
! by Jianqi Huang 2023/07/01
!=========================================================================
  use variables, only : crys_sys, spa_ord, spa_ele, &
                        natom_prim, a_prim, b_prim, iso_pg, stdout
  use mat_trans, only : mat_inv
  use constants, only : tpi
  implicit none
!-------------------------------------------------------------------------
  integer(kind=4) :: i, j, lit_ord, lit_ord_t, nirr, nirr_t
  integer(kind=4),allocatable :: lit_spa(:), lit_spa_t(:), lit_axis(:,:), lit_axis_t(:,:)
  real(kind=8) :: k_crys(3)
  real(kind=8),allocatable :: lit_ele(:,:,:), lit_ele_t(:,:,:)
  complex(kind=8),allocatable :: lit_cha(:,:), lit_cha_t(:,:)
  complex(kind=8),allocatable :: lit_irr(:,:,:,:), lit_irr_t(:,:,:,:)
  character(len=3) :: lit_cg
  character(len=2),allocatable :: lit_symb(:), lit_symb_t(:)
!-------------------------------------------------------------------------
  !
  open(unit=stdout, file='IrrepBZ.txt',status='replace')
  !
  ! read in
  call read_poscar()
  !
  ! primitive unit vectors
  call primitive_cell()
  !
  ! crystal system
  call holohedry_group()
  write (stdout,*) 'crystal system: ', crys_sys
  write (stdout,*)
  !
  ! space group elements
  call space_group()
  !
  ! isogonal point group
  call isogonal_pg(spa_ord,spa_ele,iso_pg)
  write (stdout,*) 'isogonal point group is: ', iso_pg
  write (stdout,*)
  write (stdout,*)
  !
  ! reciprocal unit vectors
  b_prim = transpose(mat_inv(a_prim))*tpi
  write (stdout,*) 'reciprocal basis vectors (Cart) :'
  write (stdout,'(3f15.10)') b_prim(:,1)
  write (stdout,'(3f15.10)') b_prim(:,2)
  write (stdout,'(3f15.10)') b_prim(:,3)
  write (stdout,*)
  !
  write (stdout, '(a/)') REPEAT('=', 100)
  !
  !call qe_band()
  !
  ! little group
  k_crys = (/ 0.d0/3, 0.d0/6, 0.d0/3 /)
  write (stdout,'("k-point (crys): ",3f15.10)') k_crys
  !
  ! without TRS
  call little_order(k_crys, lit_ord)
  write (stdout,'(a,i3)') 'the order of little co-group is:', lit_ord
  allocate(lit_spa(lit_ord),lit_ele(3,4,lit_ord), lit_symb(lit_ord), lit_axis(3,lit_ord))
  call little_group(k_crys, lit_ord, lit_spa, lit_ele, lit_symb, lit_axis, lit_cg)
  write (stdout,*) 'little co-group of given k-point is: ', lit_cg
  write (stdout,*)
  do i = 1, lit_ord
    write (stdout,'(5x,a1,i2.2,a5,3i3)') '#', i, lit_symb(i), lit_axis(:,i)
    write (stdout,'(3f15.10,f20.10)') lit_ele(1,1:3,i), lit_ele(1,4,i)
    write (stdout,'(3f15.10,f20.10)') lit_ele(2,1:3,i), lit_ele(2,4,i)
    write (stdout,'(3f15.10,f20.10)') lit_ele(3,1:3,i), lit_ele(3,4,i)
    write (stdout,*)
  end do
  write (stdout,*)
  !
  ! small vector representation
  write (stdout,*) 'small vector representations:'
  allocate(lit_irr(nint(sqrt(dble(lit_ord))),nint(sqrt(dble(lit_ord))),lit_ord,lit_ord))
  allocate(lit_cha(lit_ord,lit_ord))
  call small_rep(k_crys,lit_ord,lit_ele,lit_symb,lit_axis,nirr,lit_irr,lit_cha,.false.)
  !
  write (stdout,'(a,i3)') 'number of irreps:', nirr
  do i = 1, lit_ord
    write (stdout,'(7x,a,8x)',advance='no') lit_symb(i)
  end do
  write (stdout,*)
  do i = 1, nirr
    write (stdout,'(5x,a,i2,5x,a,i2)') '#',i,'Spacial deg:', nint(real(lit_cha(1,i)))
    do j = 1, lit_ord
      write (stdout,'(2x,"(",f6.3,","f6.3,")",$)') lit_cha(j,i)
    end do
    write (stdout,*)
    write (stdout,*)
  end do
  !
  deallocate(lit_spa, lit_symb, lit_axis, lit_ele, lit_irr, lit_cha)
  write (stdout,*)
  !
  ! TRS
  call little_order_t(k_crys, lit_ord, lit_ord_t)
  write (stdout,'(a,i3)') 'the order of little co-group is:', lit_ord_t
  allocate(lit_spa_t(lit_ord_t),lit_ele_t(3,4,lit_ord_t),lit_symb_t(lit_ord_t),lit_axis_t(3,lit_ord_t))
  call little_group_t(k_crys, lit_ord, lit_ord_t, lit_spa_t, lit_ele_t, lit_symb_t, lit_axis_t)
  allocate(lit_irr_t(nint(sqrt(dble(lit_ord_t))),nint(sqrt(dble(lit_ord_t))),lit_ord_t,lit_ord_t))
  allocate(lit_cha_t(lit_ord_t,lit_ord_t))
  call small_rep_t(k_crys,lit_ord,lit_ord_t,lit_ele_t,lit_symb_t,lit_axis_t,nirr_t,lit_irr_t,lit_cha_t,.false.)
  write (stdout,'(a,i3)') 'number of irreps:', nirr_t
  do i = 1, lit_ord_t
    write (stdout,'(7x,a,8x)',advance='no') lit_symb_t(i)
  end do
  write (stdout,*)
  do i = 1, nirr_t
    write (stdout,'(5x,a,i2,5x,a,i2)') '#',i,'Spacial deg:', nint(real(lit_cha_t(1,i)))
    do j = 1, lit_ord_t
      write (stdout,'(2x,"(",f6.3,","f6.3,")",$)') lit_cha_t(j,i)
    end do
    write (stdout,*)
    write (stdout,*)
  end do
  !
  deallocate(lit_spa_t, lit_symb_t, lit_axis_t, lit_ele_t, lit_irr_t, lit_cha_t)
  write (stdout,*)
  !
  write (stdout, '(a/)') REPEAT('=', 100)
  !
    ! SOC~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ! without TRS
  call little_order(k_crys, lit_ord)
  write (stdout,'(a,i3)') 'the order of little co-group is:', lit_ord
  allocate(lit_spa(lit_ord),lit_ele(3,4,lit_ord), lit_symb(lit_ord), lit_axis(3,lit_ord))
  call little_group(k_crys, lit_ord, lit_spa, lit_ele, lit_symb, lit_axis, lit_cg)
  write (stdout,*) 'little co-group of given k-point is: ', lit_cg
  write (stdout,*)
  !
  ! small spinor representation
  write (stdout,*) 'small vector representations:'
  allocate(lit_irr(nint(sqrt(dble(lit_ord))),nint(sqrt(dble(lit_ord))),lit_ord,lit_ord))
  allocate(lit_cha(lit_ord,lit_ord))
  call small_rep(k_crys,lit_ord,lit_ele,lit_symb,lit_axis,nirr,lit_irr,lit_cha,.true.)
  !
  write (stdout,'(a,i3)') 'number of irreps:', nirr
  do i = 1, lit_ord
    write (stdout,'(7x,a,8x)',advance='no') lit_symb(i)
  end do
  write (stdout,*)
  do i = 1, nirr
    write (stdout,'(5x,a,i2,5x,a,i2)') '#',i,'Spacial deg:', nint(real(lit_cha(1,i)))
    do j = 1, lit_ord
      write (stdout,'(2x,"(",f6.3,","f6.3,")",$)') lit_cha(j,i)
    end do
    write (stdout,*)
    write (stdout,*)
  end do
  !
  deallocate(lit_spa, lit_symb, lit_axis, lit_ele, lit_irr, lit_cha)
  write (stdout,*)
  !
  ! TRS
  call little_order_t(k_crys, lit_ord, lit_ord_t)
  write (stdout,'(a,i3)') 'the order of little co-group is:', lit_ord_t
  allocate(lit_spa_t(lit_ord_t),lit_ele_t(3,4,lit_ord_t))
  allocate(lit_symb_t(lit_ord_t))
  allocate(lit_axis_t(3,lit_ord_t))
  call little_group_t(k_crys, lit_ord, lit_ord_t, lit_spa_t, lit_ele_t, lit_symb_t, lit_axis_t)
  allocate(lit_irr_t(nint(sqrt(dble(lit_ord_t))),nint(sqrt(dble(lit_ord_t))),lit_ord_t,lit_ord_t))
  allocate(lit_cha_t(lit_ord_t,lit_ord_t))
  call small_rep_t(k_crys,lit_ord,lit_ord_t,lit_ele_t,lit_symb_t,lit_axis_t,nirr_t,lit_irr_t,lit_cha_t,.true.)
  write (stdout,'(a,i3)') 'number of irreps:', nirr_t
  do i = 1, lit_ord_t
    write (stdout,'(7x,a,8x)',advance='no') lit_symb_t(i)
  end do
  write (stdout,*)
  do i = 1, nirr_t
    write (stdout,'(5x,a,i2,5x,a,i2)') '#',i,'Spacial deg:', nint(real(lit_cha_t(1,i)))
    do j = 1, lit_ord_t
      write (stdout,'(2x,"(",f6.3,","f6.3,")",$)') lit_cha_t(j,i)
    end do
    write (stdout,*)
    write (stdout,*)
  end do
  !
  deallocate(lit_spa_t, lit_symb_t, lit_axis_t, lit_ele_t, lit_irr_t, lit_cha_t)
  write (stdout,*)
  !
  close(stdout)
end program