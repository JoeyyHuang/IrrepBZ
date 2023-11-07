subroutine read_poscar()
!=========================================================================
! read the crystal structure in POSCAR file with atomic positions listed in
! crystal or Cartesian frame
!=========================================================================
  use variables, only : ntype, natom_conv_type, symb_type, a_conv, natom_conv, &
                        atom_conv_crys, atom_conv_cart, omega_conv, stdout
  use mat_trans, only : mat_inv
  implicit none
!-------------------------------------------------------------------------
! local variables  
  integer(kind=4) :: i, p
  real(kind=8) :: fac_lat ! lattice factor
  character(len=100) type_str, frame
!------------------------------------------------------------------------- 
  !
  open(unit=10, file='POSCAR', status='old')
  read(10,*)
  read(10,*) fac_lat
  do i = 1, 3
    read(10,*) a_conv(:,i)
  end do
  a_conv = fac_lat * a_conv
  !
  ! get the ntype
  p = 1
  read(10,'(a100)') type_str
  do while(.true.)
    i = verify(type_str(p:), ' ')  !-- Find next non-blank
    if (i == 0) exit               !-- No word found
    ntype = ntype + 1              !-- Found something
    p = p + i - 1                  !-- Move to start of the word
    i = scan(type_str(p:), ' ')    !-- Find next blank
    if (i == 0) exit               !-- No blank found
    p = p + i - 1                  !-- Move to the blank
  end do
  !
  backspace(10)
  allocate(symb_type(ntype),natom_conv_type(ntype))
  read(10,*) symb_type
  read(10,*) natom_conv_type
  natom_conv = sum(natom_conv_type)
  allocate(atom_conv_crys(3,natom_conv),atom_conv_cart(3,natom_conv))
  read(10,'(a100)') frame
  if(trim(adjustl(frame)) == 'D' .or. trim(adjustl(frame)) == 'd') then
    do i = 1, natom_conv
      read(10,*) atom_conv_crys(:,i)
      atom_conv_cart(:,i) = matmul(a_conv,atom_conv_crys(:,i))
    end do
  else if(trim(adjustl(frame)) == 'C' .or. trim(adjustl(frame)) == 'c') then
    do i = 1, natom_conv
      read(10,*) atom_conv_cart(:,i)
      atom_conv_crys(:,i) = matmul(mat_inv(a_conv),atom_conv_cart(:,i))
    end do
  else
    write (stdout,*) 'wrong frame option line in POSCAR, "D" or "d" for crystal frame and "C" or "c" for Cartesian'
    stop
  end if
  close(10)
  !
  omega_conv = (a_conv(2,1)*a_conv(3,2)-a_conv(2,2)*a_conv(3,1))*a_conv(1,3) - &
               (a_conv(1,1)*a_conv(3,2)-a_conv(1,2)*a_conv(3,1))*a_conv(2,3) + &
               (a_conv(1,1)*a_conv(2,2)-a_conv(1,2)*a_conv(2,1))*a_conv(3,3)
  !
  return
end subroutine read_poscar