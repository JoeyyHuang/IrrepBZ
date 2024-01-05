# Makefile of IrrepBZ code
# replace FC with your appropriate Fortran compiler if needed

FC = gfortran
obj = constants.o Herring.o holohedry_group.o IrrepBZ.o isogonal_pg.o \
      little_group.o little_group_TRS.o matrix_trans.o primitive_cell.o \
      qe_bands.o qe_phonons.o QR.o read_binary.o read_poscar.o read_scf.o \
      rot_axis.o small_rep.o small_rep_TRS.o space_group.o spinor.o spinor_TRS.o \
      SU2.o variables.o vasp_bands.o vasp_phonons.o

IrrepBZ: $(obj)
	$(FC) -o IrrepBZ $(obj)
constants.o constants.mod: ./src/constants.f90
	$(FC) -c ./src/constants.f90
Herring.o: ./src/Herring.f90 constants.mod mat_trans.mod variables.mod
	$(FC) -c ./src/Herring.f90
holohedry_group.o: ./src/holohedry_group.f90 constants.mod mat_trans.mod variables.mod
	$(FC) -c ./src/holohedry_group.f90
IrrepBZ.o: ./src/IrrepBZ.f90 constants.mod mat_trans.mod variables.mod
	$(FC) -c ./src/IrrepBZ.f90
isogonal_pg.o: ./src/isogonal_pg.f90 constants.mod mat_trans.mod
	$(FC) -c ./src/isogonal_pg.f90
little_group.o: ./src/little_group.f90 constants.mod mat_trans.mod variables.mod
	$(FC) -c ./src/little_group.f90
little_group_TRS.o: ./src/little_group_TRS.f90 constants.mod mat_trans.mod variables.mod
	$(FC) -c ./src/little_group_TRS.f90
matrix_trans.o mat_trans.mod: ./src/matrix_trans.f90
	$(FC) -c ./src/matrix_trans.f90
primitive_cell.o: ./src/primitive_cell.f90 constants.mod mat_trans.mod variables.mod
	$(FC) -c ./src/primitive_cell.f90
qe_bands.o: ./src/qe_bands.f90 constants.mod mat_trans.mod read_binary.mod variables.mod
	$(FC) -c ./src/qe_bands.f90
qe_phonons.o: ./src/qe_phonons.f90 constants.mod mat_trans.mod read_binary.mod variables.mod
	$(FC) -c ./src/qe_phonons.f90
QR.o: ./src/QR.f90
	$(FC) -c ./src/QR.f90
read_binary.o read_binary.mod: ./src/read_binary.f90
	$(FC) -c ./src/read_binary.f90
read_poscar.o: ./src/read_poscar.f90 constants.mod mat_trans.mod variables.mod
	$(FC) -c ./src/read_poscar.f90
read_scf.o: ./src/read_scf.f90 constants.mod mat_trans.mod variables.mod
	$(FC) -c ./src/read_scf.f90
rot_axis.o: ./src/rot_axis.f90 constants.mod mat_trans.mod variables.mod
	$(FC) -c ./src/rot_axis.f90
small_rep.o: ./src/small_rep.f90 constants.mod
	$(FC) -c ./src/small_rep.f90
small_rep_TRS.o: ./src/small_rep_TRS.f90 constants.mod
	$(FC) -c ./src/small_rep_TRS.f90
space_group.o: ./src/space_group.f90 constants.mod mat_trans.mod variables.mod
	$(FC) -c ./src/space_group.f90
spinor.o: ./src/spinor.f90 constants.mod variables.mod
	$(FC) -c ./src/spinor.f90
spinor_TRS.o: ./src/spinor_TRS.f90 constants.mod variables.mod
	$(FC) -c ./src/spinor_TRS.f90
SU2.o: ./src/SU2.f90 constants.mod mat_trans.mod variables.mod
	$(FC) -c ./src/SU2.f90
variables.o variables.mod: ./src/variables.f90
	$(FC) -c ./src/variables.f90
vasp_bands.o: ./src/vasp_bands.f90 constants.mod mat_trans.mod variables.mod
	$(FC) -c ./src/vasp_bands.f90
vasp_phonons.o: ./src/vasp_phonons.f90 constants.mod mat_trans.mod variables.mod
	$(FC) -c ./src/vasp_phonons.f90

clean:
	rm -f IrrepBZ $(obj) constants.mod mat_trans.mod variables.mod read_binary.mod
