ALL: blib exec

#compilation and various flags
DIRS    = 
EXEC    = matgen.exe
CFLAGS	= 
FFLAGS	= 
CPPFLAGS	= 
FPPFLAGS	=
CLEANFILES  = ${EXEC}
OFILES= ${wildcard ./*.o}
NBLOCK = 2
DEBUG_VALGRIND = valgrind

include ${PETSC_DIR}/lib/petsc/conf/variables
include ${PETSC_DIR}/lib/petsc/conf/rules

blib :
	-@echo "BEGINNING TO COMPILE LIBRARIES "
	-@echo "========================================="
	-@${OMAKE}  PETSC_ARCH=${PETSC_ARCH} PETSC_DIR=${PETSC_DIR} ACTION=libfast tree
	-@echo "Completed building libraries"
	-@echo "========================================="

distclean :
	-@echo "Cleaning application and libraries"
	-@echo "========================================="
	-@${OMAKE} PETSC_ARCH=${PETSC_ARCH}  PETSC_DIR=${PETSC_DIR} clean
	-${RM} ${OFILES}
	-@echo "Finised cleaning application and libraries"
	-@echo "========================================="	

exec: edmg.o ebmg.o libs.o main.o
	-@echo "BEGINNING TO COMPILE APPLICATION "
	-@echo "========================================="
	-@${CLINKER} -o ${EXEC} edmg.o ebmg.o libs.o main.o ${PETSC_LIB}
	-@echo "Completed building application"
	-@echo "========================================="


run1:
	-@${MPIEXEC} -np 3 ./matgen.exe -vfile lsqr.bin -n 274 -nzeros 100 -dense
run2:
	-@${MPIEXEC} -np 40 ./matgen.exe -n 1000000 -nzeros 999999
run3:
	-@${MPIEXEC} -np 2 ./matgen.exe -n 300 -nzeros 200 -realMat -edmg
run4:
	-@${MPIEXEC} -np 1 ./matgen.exe -vfile lsqr.bin -n 274 -nzeros 200 -realMat

run5:
	./matgen.exe -n 330001 -nzeros 330000

run6:
	-@${MPIEXEC} -np 2 ./matgen.exe -n 274 -vfile lsqr.bin
