subroutine read_scf()
!=========================================================================
! read the crystal structure in scf.in file with atomic positions in
! crystal or Cartesian frame
!=========================================================================
  use variables, only : ntype, natom_conv_type, symb_type, a_conv, natom_conv, &
                        atom_conv_crys, atom_conv_cart, omega_conv
  use mat_trans, only : mat_inv
  use constants, only : bohr2ang
  implicit none
!-------------------------------------------------------------------------
! local variables  
  integer(kind=4) :: i, i1, j1, i2, j2, i3, j3, i4, j4, p, ibrav, nat, ntyp
  real(kind=8) :: celldm(6), a1(3), a2(3), a3(3)
  character(len=132) :: line
  character(len=4),allocatable :: symb(:)
  namelist /SYSTEM/ ibrav, celldm, nat, ntyp
!------------------------------------------------------------------------- 
  ibrav = 0
  celldm = 0.d0
  a1 = 0.d0; a2 = 0.d0; a3 = 0.d0
  open(unit=10, file='scf.in', status='old')
  read (10,nml=SYSTEM)
  close(10)
  celldm(1) = celldm(1)*bohr2ang  ! convert to angstrom units
  !
  if(ibrav == 0) then
    open(unit=10, file='scf.in', status='old')
    do while(.true.)
      read(10,'(a132)') line
      p = index(line,'CELL_PARAMETERS',back=.false.,kind=4)
      if(p /= 0) then
        i1 = index(line,'alat',back=.false.,kind=4)
        j1 = index(line,'ALAT',back=.false.,kind=4)
        i2 = index(line,'bohr',back=.false.,kind=4)
        j2 = index(line,'BOHR',back=.false.,kind=4)
        i3 = index(line,'angstrom',back=.false.,kind=4)
        j3 = index(line,'ANGSTROM',back=.false.,kind=4)
        if(i1 /= 0 .or. j1 /= 0) then
          read(10,*) a1; a1 = a1*celldm(1)
          read(10,*) a2; a2 = a2*celldm(1)
          read(10,*) a3; a3 = a3*celldm(1)
        elseif(i2 /= 0 .or. j2 /= 0) then
          read(10,*) a1; a1 = a1*bohr2ang
          read(10,*) a2; a2 = a2*bohr2ang
          read(10,*) a3; a3 = a3*bohr2ang
        elseif(i3 /= 0 .or. j3 /= 0) then
          read(10,*) a1
          read(10,*) a2
          read(10,*) a3
        end if
        close(10)
        exit
      end if
    end do
  end if
  !
  call bravais(ibrav, celldm, a1, a2, a3)
  !
  a_conv(:,1) = a1(:)
  a_conv(:,2) = a2(:)
  a_conv(:,3) = a3(:)
  omega_conv = (a_conv(2,1)*a_conv(3,2)-a_conv(2,2)*a_conv(3,1))*a_conv(1,3) - &
               (a_conv(1,1)*a_conv(3,2)-a_conv(1,2)*a_conv(3,1))*a_conv(2,3) + &
               (a_conv(1,1)*a_conv(2,2)-a_conv(1,2)*a_conv(2,1))*a_conv(3,3)
  !
  ntype = ntyp; natom_conv = nat
  allocate(symb(natom_conv),symb_type(ntype),natom_conv_type(ntype))
  allocate(atom_conv_crys(3,natom_conv),atom_conv_cart(3,natom_conv))
  !
  open(unit=10, file='scf.in', status='old')
  do while(.true.)
    read(10,'(a132)') line
    p = index(line,'ATOMIC_POSITIONS',back=.false.,kind=4)
    if(p /= 0) then
      i1 = index(line,'alat',back=.false.,kind=4)
      j1 = index(line,'ALAT',back=.false.,kind=4)
      i2 = index(line,'bohr',back=.false.,kind=4)
      j2 = index(line,'BOHR',back=.false.,kind=4)
      i3 = index(line,'angstrom',back=.false.,kind=4)
      j3 = index(line,'ANGSTROM',back=.false.,kind=4)
      i4 = index(line,'crystal',back=.false.,kind=4)
      j4 = index(line,'CRYSTAL',back=.false.,kind=4)
      if(i1 /= 0 .or. j1 /= 0) then
        do i = 1, natom_conv
          read(10,*) symb(i), atom_conv_cart(:,i)
          atom_conv_cart(:,i) = atom_conv_cart(:,i)*celldm(1)
          atom_conv_crys(:,i) = matmul(mat_inv(a_conv),atom_conv_cart(:,i))
        end do
      elseif(i2 /= 0 .or. j2 /= 0) then
        do i = 1, natom_conv
          read(10,*) symb(i), atom_conv_cart(:,i)
          atom_conv_cart(:,i) = atom_conv_cart(:,i)*bohr2ang
          atom_conv_crys(:,i) = matmul(mat_inv(a_conv),atom_conv_cart(:,i))
        end do
      elseif(i3 /= 0 .or. j3 /= 0) then
        do i = 1, natom_conv
          read(10,*) symb(i), atom_conv_cart(:,i)
          atom_conv_crys(:,i) = matmul(mat_inv(a_conv),atom_conv_cart(:,i))
        end do
      elseif(i4 /= 0 .or. j4 /= 0) then
        do i = 1, natom_conv
          read(10,*) symb(i), atom_conv_crys(:,i)
          atom_conv_cart(:,i) = matmul(a_conv,atom_conv_crys(:,i))
        end do
      end if
      close(10)
      exit
    end if
  end do
  !
  p = 1
  symb_type = '????'
  natom_conv_type = 0
  symb_type(p) = symb(p)
  natom_conv_type(p) = 1
  if(natom_conv > 1) then
    do i = 2, natom_conv
      if(symb(i) == symb_type(p)) then
        natom_conv_type(p) = natom_conv_type(p) + 1
      else
        p = p + 1
        symb_type(p) = symb(i)
        natom_conv_type(p) = natom_conv_type(p) + 1
      end if
    end do
  end if
  deallocate(symb)
  !
  return
end subroutine read_scf
!
SUBROUTINE bravais(ibrav, celldm, a1, a2, a3)
!=========================================================================
  !     Adopted from latgen.f90 in Modules of Quantum Espresso package
  !
  !     sets up the crystallographic vectors a1, a2, and a3.
  !
  !     ibrav is the structure index:
  !       1  cubic P (sc)                8  orthorhombic P
  !       2  cubic F (fcc)               9  1-face (C) centered orthorhombic
  !       3  cubic I (bcc)              10  all face centered orthorhombic
  !       4  hexagonal and trigonal P   11  body centered orthorhombic
  !       5  trigonal R, 3-fold axis c  12  monoclinic P (unique axis: c)
  !       6  tetragonal P (st)          13  one face (base) centered monoclinic
  !       7  tetragonal I (bct)         14  triclinic P
  !     Also accepted:
  !       0  "free" structure          -12  monoclinic P (unique axis: b)
  !      -3  cubic bcc with a more symmetric choice of axis
  !      -5  trigonal R, threefold axis along (111)
  !      -9  alternate description for base centered orthorhombic
  !     -13  one face (base) centered monoclinic (unique axis: b)
  !      91  1-face (A) centered orthorombic
  !
  !     celldm are parameters which fix the shape of the unit cell
  !
  !     NOTA BENE: all axis sets are right-handed
!=========================================================================
  use constants, only : tol, stdout
  IMPLICIT NONE
  INTEGER(kind=4), INTENT(in) :: ibrav
  real(kind=8), INTENT(inout) :: celldm(6)
  real(kind=8), INTENT(inout) :: a1(3), a2(3), a3(3)
  !
  real(kind=8), PARAMETER:: sr2 = sqrt(2.d0), sr3 = sqrt(3.d0)
  INTEGER(kind=4) :: ir
  real(kind=8) :: term, cbya, term1, term2, singam, sen
  !
  !  user-supplied lattice vectors
  IF (ibrav /= 0) THEN
     a1(:) = 0.d0
     a2(:) = 0.d0
     a3(:) = 0.d0
  ENDIF
  !
  !  index of bravais lattice supplied
  IF (ibrav == 1) THEN   ! simple cubic lattice
     a1(1)=celldm(1)
     a2(2)=celldm(1)
     a3(3)=celldm(1)
  ELSEIF (ibrav == 2) THEN   ! fcc lattice
     term=celldm(1)/2.d0
     a1(1)=-term
     a1(3)=term
     a2(2)=term
     a2(3)=term
     a3(1)=-term
     a3(2)=term
  ELSEIF (abs(ibrav) == 3) THEN   ! bcc lattice
     term=celldm(1)/2.d0
     DO ir=1,3
        a1(ir)=term
        a2(ir)=term
        a3(ir)=term
     ENDDO
     IF ( ibrav < 0 ) THEN
        a1(1)=-a1(1)
        a2(2)=-a2(2)
        a3(3)=-a3(3)
     ELSE
        a2(1)=-a2(1)
        a3(1)=-a3(1)
        a3(2)=-a3(2)
     ENDIF
  ELSEIF (ibrav == 4) THEN   ! hexagonal lattice
     cbya=celldm(3)
     a1(1)=celldm(1)
     a2(1)=-celldm(1)/2.d0
     a2(2)=celldm(1)*sr3/2.d0
     a3(3)=celldm(1)*cbya
  ELSEIF (abs(ibrav) == 5) THEN   ! trigonal lattice
     term1=sqrt(1.d0 + 2.d0*celldm(4))
     term2=sqrt(1.d0 - celldm(4))
     IF ( ibrav == 5) THEN   ! threefold axis along c (001)
        a2(2)=sr2*celldm(1)*term2/sr3
        a2(3)=celldm(1)*term1/sr3
        a1(1)=celldm(1)*term2/sr2
        a1(2)=-a1(1)/sr3
        a1(3)= a2(3)
        a3(1)=-a1(1)
        a3(2)= a1(2)
        a3(3)= a2(3)
     ELSEIF ( ibrav == -5) THEN   ! threefold axis along (111)
        ! Notice that in the cubic limit (alpha=90, celldm(4)=0, term1=term2=1)
        ! does not yield the x,y,z axis, but an equivalent rotated triplet:
        !    a/3 (-1,2,2), a/3 (2,-1,2), a/3 (2,2,-1)
        ! If you prefer the x,y,z axis as cubic limit, you should modify the
        ! definitions of a1(1) and a1(2) as follows:'
        !    a1(1) = celldm(1)*(term1+2.0_dp*term2)/3.0_dp
        !    a1(2) = celldm(1)*(term1-term2)/3.0_dp
        ! (info by G. Pizzi and A. Cepellotti)
        a1(1) = celldm(1)*(term1-2.d0*term2)/3.d0
        a1(2) = celldm(1)*(term1+term2)/3.d0
        a1(3) = a1(2)
        a2(1) = a1(3)
        a2(2) = a1(1)
        a2(3) = a1(2)
        a3(1) = a1(2)
        a3(2) = a1(3)
        a3(3) = a1(1)
     ENDIF
  ELSEIF (ibrav == 6) THEN   ! tetragonal lattice
     cbya=celldm(3)
     a1(1)=celldm(1)
     a2(2)=celldm(1)
     a3(3)=celldm(1)*cbya
  ELSEIF (ibrav == 7) THEN   ! body centered tetragonal lattice
     cbya=celldm(3)
     a2(1)=celldm(1)/2.d0
     a2(2)=a2(1)
     a2(3)=cbya*celldm(1)/2.d0
     a1(1)= a2(1)
     a1(2)=-a2(1)
     a1(3)= a2(3)
     a3(1)=-a2(1)
     a3(2)=-a2(1)
     a3(3)= a2(3)
  ELSEIF (ibrav == 8) THEN   ! Simple orthorhombic lattice
     a1(1)=celldm(1)
     a2(2)=celldm(1)*celldm(2)
     a3(3)=celldm(1)*celldm(3)
  ELSEIF ( abs(ibrav) == 9) THEN   ! One face (base) centered orthorhombic lattice  (C type)
     IF ( ibrav == 9 ) THEN   ! old PWscf description
        a1(1) = 0.5d0 * celldm(1)
        a1(2) = a1(1) * celldm(2)
        a2(1) = - a1(1)
        a2(2) = a1(2)
     ELSE   ! alternate description
        a1(1) = 0.5d0 * celldm(1)
        a1(2) =-a1(1) * celldm(2)
        a2(1) = a1(1)
        a2(2) =-a1(2)
     ENDIF
     a3(3) = celldm(1) * celldm(3)
  ELSEIF ( ibrav == 91 ) THEN   ! One face (base) centered orthorhombic lattice  (A type)
     a1(1) = celldm(1)
     a2(2) = celldm(1) * celldm(2) * 0.5d0
     a2(3) = - celldm(1) * celldm(3) * 0.5d0
     a3(2) = a2(2)
     a3(3) = - a2(3)
  ELSEIF (ibrav == 10) THEN   ! All face centered orthorhombic lattice
     a2(1) = 0.5d0 * celldm(1)
     a2(2) = a2(1) * celldm(2)
     a1(1) = a2(1)
     a1(3) = a2(1) * celldm(3)
     a3(2) = a2(1) * celldm(2)
     a3(3) = a1(3)
  ELSEIF (ibrav == 11) THEN   ! Body centered orthorhombic lattice
     a1(1) = 0.5d0 * celldm(1)
     a1(2) = a1(1) * celldm(2)
     a1(3) = a1(1) * celldm(3)
     a2(1) = - a1(1)
     a2(2) = a1(2)
     a2(3) = a1(3)
     a3(1) = - a1(1)
     a3(2) = - a1(2)
     a3(3) = a1(3)
  ELSEIF (ibrav == 12) THEN   ! Simple monoclinic lattice, unique (i.e. orthogonal to a) axis: c
     sen=sqrt(1.d0-celldm(4)**2)
     a1(1)=celldm(1)
     a2(1)=celldm(1)*celldm(2)*celldm(4)
     a2(2)=celldm(1)*celldm(2)*sen
     a3(3)=celldm(1)*celldm(3)
  ELSEIF (ibrav ==-12) THEN   ! Simple monoclinic lattice, unique axis: b (more common)
     sen=sqrt(1.d0-celldm(5)**2)
     a1(1)=celldm(1)
     a2(2)=celldm(1)*celldm(2)
     a3(1)=celldm(1)*celldm(3)*celldm(5)
     a3(3)=celldm(1)*celldm(3)*sen
  ELSEIF (ibrav == 13) THEN   ! One face centered monoclinic lattice unique axis c
     sen = sqrt( 1.d0 - celldm(4) ** 2 )
     a1(1) = 0.5d0 * celldm(1)
     a1(3) =-a1(1) * celldm(3)
     a2(1) = celldm(1) * celldm(2) * celldm(4)
     a2(2) = celldm(1) * celldm(2) * sen
     a3(1) = a1(1)
     a3(3) =-a1(3)
  ELSEIF (ibrav == -13) THEN   ! One face centered monoclinic lattice unique axis b
     sen = sqrt( 1.d0 - celldm(5) ** 2 )
     a1(1) = 0.5d0 * celldm(1)
     a1(2) =-a1(1) * celldm(2)
     a2(1) = a1(1)
     a2(2) =-a1(2)
     a3(1) = celldm(1) * celldm(3) * celldm(5)
     a3(3) = celldm(1) * celldm(3) * sen
  ELSEIF (ibrav == 14) THEN   ! Triclinic lattice
     singam=sqrt(1.d0-celldm(6)**2)
     term= (1.d0+2.d0*celldm(4)*celldm(5)*celldm(6)-celldm(4)**2-celldm(5)**2-celldm(6)**2)
     IF (term < tol) write (stdout,'(a)') 'celldm do not make sense, check your data in scf.in'
     term= sqrt(term/(1.d0-celldm(6)**2))
     a1(1)=celldm(1)
     a2(1)=celldm(1)*celldm(2)*celldm(6)
     a2(2)=celldm(1)*celldm(2)*singam
     a3(1)=celldm(1)*celldm(3)*celldm(5)
     a3(2)=celldm(1)*celldm(3)*(celldm(4)-celldm(5)*celldm(6))/singam
     a3(3)=celldm(1)*celldm(3)*term
  ELSE
     write (stdout,'(a)') 'Nonexistent bravais lattice in scf.in'
  ENDIF
  !
  RETURN
END SUBROUTINE bravais