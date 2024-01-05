program IrrepBZ
!=========================================================================
! Functionality: 
! 1. If the crystal structure is only present either in the POSCAR file (VASP) or scf.in file (Quantum Espresso), 
!    the information related to spatial symmetry such as crystal system, holohedral point group, 
!    symmetry operations of the space group, and so on, will be documented in the 'IrrepBZ.txt' file.
!
! 2. besides crystal structure file, if a 'k-points-list' file is also present, which provides the coordinates of the 
!    k-points to be investigated, then the small representaions of k-point little groups, under different scenarios, 
!    considering the presence or absence of spin-orbit coupling and time-reversal symmetry, will be generated on-the-fly 
!    in seperate 'Irrep_ik*.txt' files, one for each k-point.
!
! 3. besides crystal structure file, if there is a 'qe_band' folder, consisting of 'wfc*.dat' files along with 'data-file-schema.xml' file 
!    obtained from the prefix.save folder after performing a nscf calculation with Quantum Espresso, the irreducible representations of 
!    electronic bands present in, like, wfc1.dat will be saved to a 'band_ik1.txt' file within the same folder.
!
! 4. besides crystal structure file, if there is an additional 'qe_phonon' folder, containing 'matdyn.freq' and 'matdyn.eig' files 
!    generated through the implementation of matdyn.x in the phonon calculation utilizing Quantum Espresso, the irreducible representations
!    of normal vibrations at specific N-th q-point will be recorded in a 'phonon_iqN.txt' file within the same folder.
!
! 5. besides crystal structure file, if there is also a 'vasp_band' folder which includes the 'WAVECAR' file obtained from a non-consistent 
!    calculation using VASP, the irreducible representations of electronic bands present in the WAVECAR file will be saved in seperate 
!    'band_ik*.txt' files within the same folder, each corresponding to one k-point.
!
! 6. besides crystal structure file, if there is also a 'vasp_phonon' folder that includes the 'band.yaml' or 'qpoints.yaml' files 
!    obtained from a phonon calculation using the Phonopy package interfaced with VASP, the irreducible representations of normal vibrations 
!    at specific N-th q-point will be recorded in a 'phonon_iqN.txt' file within the same folder.
!
! Beware:
! * Please note that only one crystal structure is processed at a time, so either the POSCAR or scf.in file can exist, but not both simultaneously.
!
! * In the scf.in file, only the tags of ibrav, celldm, nat, and ntyp are allowed to be present in the SYSTEM namelist.
!
! * If the scf.in file exists, it is optional for the 'k-points-list' file, 'qe_band' folder, and 'qe_phonon' folder to exist. 
!   They can all exist together, none of them can exist, or they can exist in any combination. 
!   The same rule applies to the 'k-points-list' file, 'vasp_band' folder, and 'vasp_phonon' folder in the case of the POSCAR file.
!
! * Regarding the 'vasp_phonon' folder, the q-points can be specified either in a band format or as discrete q-points in the Phonopy implementation. 
!   These correspond respectively to the resulting files of 'band.yaml' and 'qpoints.yaml'. 
!   Therefore, only one of these files is required and they cannot coexist together.
!
! 
! by Jianqi Huang, July 2023
!=========================================================================
  use variables, only : natom_conv_type,atom_conv_crys,atom_conv_cart,symb_type,natom_prim_type,&
                        atom_prim_crys,atom_prim_cart,holo_ele,spa_ele,spa_axis,spa_symb,b_prim
  use mat_trans, only : mat_inv
  use constants, only : stdout
  implicit none
!-------------------------------------------------------------------------
  integer(kind=4) :: ik, nk, lit_ord, lit_ord_t, nirr, nirr_t
  integer(kind=4),allocatable :: lit_spa(:), lit_spa_t(:), lit_axis(:,:), lit_axis_t(:,:)
  real(kind=8) :: k_crys(3), k_cart(3)
  real(kind=8),allocatable :: lit_ele(:,:,:), lit_ele_t(:,:,:)
  complex(kind=8),allocatable :: lit_cha(:,:), lit_cha_t(:,:)
  character(len=2),allocatable :: lit_symb(:), lit_symb_t(:)
  character(len=10) klabel
  character(len=20) fname
  character(len=100) :: frame
  logical :: exst, exst_vasp, exst_qe
!-------------------------------------------------------------------------
  open(unit=stdout, file='IrrepBZ.txt',status='replace')
  inquire(file = 'POSCAR', exist = exst_vasp)
  inquire(file = 'scf.in', exist = exst_qe)
  if(exst_vasp .and. exst_qe) then
    write (stdout,'(a)') 'Error: detect the presence of both POSCAR and scf.in,'
    write (stdout,'(a)') '       make sure only one of them exists!'
    stop
  end if
  if(exst_vasp) call read_poscar()  ! read atmoic structure in POSCAR
  if(exst_qe) call read_scf()       ! read atmoic structure in scf.in
  !
  call primitive_cell()   ! primitive unit vectors
  call holohedry_group()  ! crystal system
  call space_group()      ! symmetry operations of space group
  call SU2()              ! SU2 representation of symmetry operations
  close(stdout)
  !
  ! calculate small representation for little group at given k-points
  inquire(file = 'k-points-list', exist = exst)
  if(exst) then
    open(unit=10, file='k-points-list', status='old')
    read(10,*) nk
    read(10,*) frame
    do ik = 1, nk
      if(ik == 1) then
        write(klabel,'(a)') '1st'
        write(fname,'(a10,i1.1,a4)') 'IrrepBZ_ik', ik, '.txt'
      else if(ik == 2) then
        write(klabel,'(a)') '2nd'
        write(fname,'(a10,i1.1,a4)') 'IrrepBZ_ik', ik, '.txt'
      else if(ik == 3) then
        write(klabel,'(a)') '3rd'
        write(fname,'(a10,i1.1,a4)') 'IrrepBZ_ik', ik, '.txt'
      else if(ik < 10) then
        write(klabel,'(i1,a)') ik, 'th'
        write(fname,'(a10,i1.1,a4)') 'IrrepBZ_ik', ik, '.txt'
      else if(ik < 100) then
        write(klabel,'(i2,a)') ik, 'th'
        write(fname,'(a10,i2.2,a4)') 'IrrepBZ_ik', ik, '.txt'
      else if(ik < 1000) then
        write(klabel,'(i3,a)') ik, 'th'
        write(fname,'(a10,i3.3,a4)') 'IrrepBZ_ik', ik, '.txt'
      end if
      open(unit=stdout, file=trim(adjustl(fname)), status='replace')
      write (stdout,'(a,i3)') 'Number of k points:', nk
      if(trim(adjustl(frame)) == 'D' .or. trim(adjustl(frame)) == 'd' .or. trim(adjustl(frame)) == 'Direct') then
        read(10,*) k_crys
        write (stdout,'(a,a,a,3f10.5)') 'This is the ', trim(adjustl(klabel)), ' k-point (cryst. coord.):', k_crys
      end if
      if(trim(adjustl(frame)) == 'C' .or. trim(adjustl(frame)) == 'c' .or. trim(adjustl(frame)) == 'Cartesian') then
        read(10,*) k_cart
        write (stdout,'(a,a,a,3f10.5)') 'This is the ', trim(adjustl(klabel)), ' k-point (cart. coord.):', k_cart
        k_crys = matmul(mat_inv(b_prim),k_cart)
      end if
      !!!!!!!!!!!!!!!!!!!!! vector representation!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! ordinary little group with all unitary operations
      call little_order(k_crys, lit_ord)
      allocate(lit_spa(lit_ord),lit_ele(3,4,lit_ord), lit_symb(lit_ord), lit_axis(3,lit_ord))
      call little_group(stdout, k_crys, lit_ord, lit_spa, lit_ele, lit_symb, lit_axis)
      write (stdout, '(a/)') repeat('=', 100)
      allocate(lit_cha(lit_ord,lit_ord))
      call small_rep(stdout,k_crys,lit_ord,lit_ele,lit_symb,lit_axis,nirr,lit_cha,.false.)
      call Herring(stdout, k_crys, lit_ord, lit_spa, lit_ele, lit_symb, nirr, lit_cha, .false.)
      ! extended little group when consider the anti-unitary operation of TRS as en element
      call little_order_t(k_crys, lit_ord, lit_ord_t)
      allocate(lit_spa_t(lit_ord_t),lit_ele_t(3,4,lit_ord_t),lit_symb_t(lit_ord_t),lit_axis_t(3,lit_ord_t))
      call little_group_t(stdout, k_crys, lit_ord, lit_ord_t, lit_spa_t, lit_ele_t, lit_symb_t, lit_axis_t)
      allocate(lit_cha_t(lit_ord_t,lit_ord_t))
      call small_rep_t(stdout,k_crys,lit_ord,lit_ord_t,lit_spa_t,lit_ele_t,lit_symb_t,lit_axis_t,nirr_t,lit_cha_t,.false.)
      write (stdout, '(a/)') repeat('=', 100)
      !!!!!!!!!!!!!!!!!!!!!!! spinor representation!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! ordinary little group with all unitary operations
      call small_rep(stdout,k_crys,lit_ord,lit_ele,lit_symb,lit_axis,nirr,lit_cha,.true.)
      call Herring(stdout, k_crys, lit_ord, lit_spa, lit_ele, lit_symb, nirr, lit_cha, .true.)
      deallocate(lit_spa, lit_symb, lit_axis, lit_ele, lit_cha)
      ! extended little group when consider the anti-unitary operation of TRS as en element
      call small_rep_t(stdout,k_crys,lit_ord,lit_ord_t,lit_spa_t,lit_ele_t,lit_symb_t,lit_axis_t,nirr_t,lit_cha_t,.true.)
      write (stdout, '(a/)') repeat('=', 100)
      deallocate(lit_spa_t, lit_symb_t, lit_axis_t, lit_ele_t, lit_cha_t)
      close(stdout)
    end do
    close(10)
  end if
  !
  ! group analysis of QE electronic bands
  inquire(file = 'qe_band/data-file-schema.xml', exist = exst)
  if(exst) call qe_bands()
  !
  ! group analysis of QE phonons
  inquire(file = 'qe_phonon/matdyn.eig', exist = exst)
  if(exst) call qe_phonons()
  !
  ! group analysis of VASP electronic bands
  inquire(file = 'vasp_band/WAVECAR', exist = exst)
  if(exst) call vasp_bands()
  !
  ! group analysis of VASP phonons
  inquire(file = 'vasp_phonon/band.yaml', exist = exst)
  if(exst) call vasp_phonons()
  inquire(file = 'vasp_phonon/qpoints.yaml', exist = exst)
  if(exst) call vasp_phonons()
  !
  deallocate(natom_conv_type,atom_conv_crys,atom_conv_cart,symb_type,natom_prim_type,atom_prim_crys,atom_prim_cart, &
             holo_ele,spa_ele,spa_axis,spa_symb)
  !
end program