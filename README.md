# IrrepBZ
IrrepBZ program on-the-fly generates the irreducible representations (Irreps) of wavevector group, also known as the small representations of little group, at any given points of the Brillouin zone (BZ).
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
