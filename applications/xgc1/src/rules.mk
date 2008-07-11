$(CMD): $(OBJ) $(BINPACK_LIB)
	$(CMP) -o $(CMD) $(OPT) $(OBJ) $(LIB) -L${PETSC_DIR}/lib/${PETSC_ARCH} -lcraypetsc ${PETSC_FORTRAN_LIB} ${PETSC_KSP_LIB}

$(OBJ): module.F90

.F90.o : config.xml
	python gwrite.py config.xml temp
	$(CMP) $(OPT) -I${PETSC_DIR} -I${PETSC_DIR}/bmake/${PETSC_ARCH} -I${PETSC_DIR}/include -I${PETSC_DIR}/include/finclude ${MPI_INCLUDE} -c temp/$<

clean:
	rm -f $(CMD) *.o *.mod 


