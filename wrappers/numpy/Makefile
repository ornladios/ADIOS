MPI=n
CYTHON=n

OBJ = adios.so

ifeq ($(CYTHON), y)
	OBJ := adios.cpp $(OBJ)
endif

ifeq ($(MPI),y)
	OBJ += adios_mpi.so
	ifeq ($(CYTHON), y)
		OBJ := adios_mpi.cpp $(OBJ)
	endif
endif


default:
	@echo "ADIOS wrapper for Python/Numpy"
	@echo "Have 'adios_config' and 'python' in the path and run"
	@echo "> make python"
	@echo ""
	@echo "MPI-enabled ADIOS wrapper can be built (MPI4Py is required):"
	@echo "> make MPI=y python"
	@echo ""

python: $(OBJ)

adios.cpp:  adios.pyx
	cython --cplus adios.pyx

adios.so:  
	python setup.py build_ext -lrt

adios_mpi.cpp:  adios_mpi.pyx
	cython --cplus adios_mpi.pyx

adios_mpi.so:  
	python setup_mpi.py build

clean:
	rm -rf build

install:
	@echo "Run (you may need an admin previlege):"
	@echo "> python setup.py install"
	@echo ""
	@echo "You may specify PREFIX:"
	@echo "> python setup.py install --prefix=/dir/to/install"
	@echo ""
	@echo "You can also install in your local directory, run"
	@echo "> python setup.py install --user"

