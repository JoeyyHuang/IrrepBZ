module variables
  implicit none
  save
  ! conventional cell
  integer(kind=4) :: ntype = 0   ! number of atomic types
  integer(kind=4) :: natom_conv  ! number of atoms
  integer(kind=4),allocatable :: natom_conv_type(:) ! number of atoms each type
  real(kind=8) :: omega_conv     ! volume
  real(kind=8) :: a_conv(3,3)    ! unit vectors (coords,counter)
  real(kind=8),allocatable :: atom_conv_crys(:,:)   ! crystal atomic positions (coords,counter)
  real(kind=8),allocatable :: atom_conv_cart(:,:)   ! Cartesian atomic positions (coords,counter)
  character(len=5),allocatable :: symb_type(:)      ! chemical symbol of each type
  !
  ! primitive cell
  integer(kind=4) :: natom_prim = 0 ! number of atoms
  integer(kind=4),allocatable :: natom_prim_type(:) ! number of atoms each type
  real(kind=8) :: omega_prim        ! volume
  real(kind=8) :: a_prim(3,3)       ! unit vectors (coords,counter)
  real(kind=8),allocatable :: atom_prim_crys(:,:)   ! crystal atomic positions (coords,counter)
  real(kind=8),allocatable :: atom_prim_cart(:,:)   ! Cartesian atomic positions (coords,counter)
  !
  ! reciprocal cell
  real(kind=8) :: b_prim(3,3)  ! (coords,counter)
  !
  ! symmetry group
  integer(kind=4) :: holo_ord = 0   ! order of holohedry group
  integer(kind=4) :: spa_ord = 0    ! order of space group
  real(kind=8),allocatable :: holo_ele(:,:,:)  ! holohedry elements (3,3,n)(rotation,counter) in primitive basis
  real(kind=8),allocatable :: spa_ele(:,:,:)   ! space elements (3,4,n)(rotation,fraction,counter) in primitive basis
  integer(kind=4),allocatable :: spa_axis(:,:)   ! (vector,counter_index)rotational axes (in crys)
  character(len=2),allocatable :: spa_symb(:)   ! type of space rotations
  character(len=4) :: crys_sys = 'Tric'   ! crystal system
  character(len=3) :: iso_pg   ! isogonal point group
  !
  ! standard output
  integer(kind=4) :: stdout = 6
  
  
end module variables
  