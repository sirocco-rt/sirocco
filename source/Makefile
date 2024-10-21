# This is the makefile for the sirocco related programs
#
# usage      make [CC=...] [D] sirocco
# 			 make [CC=...] [D] all to remake all of the routines, and assure they
# 			 					are properly indented.
#
# Adding D causes the routine to be run in a way that profiling and debuggers  can be used more easily.
# Otherwise the run will be optimized to run as fast as possible. CC is an option to choose
# a different compiler other than mpicc.
#

# MATOM_VER is a temporary switch that allows one to compile the code so that it will
# use the gsl routines wheen computing alpha_sp.  It should be invoked from the commad line, e.g
# make MATOM_VER=2 sirocco (or all); the Makefile should not be altered unless we change our mind
# about the integration approach to use. The default is set to use the newer brute force technique

MATOM_VER := 1

#MPICC is now default compiler
CC = mpicc
# CC = gcc	# can use GCC either from command line or by uncommenting this
FC = g77
# FC = gfortran
NVCC = 
# NVCC = nvcc  # default (and only?) CUDA compiler
# specify any extra compiler flags here
EXTRA_FLAGS = -Wno-deprecated-non-prototype -DMATOM_VER=$(MATOM_VER)
LDFLAGS =

# Use the next line if you want to change the c compiler that mpicc uses, or
# better, set this outside the Makefile. Use unset to go back to the default cc
# export OMPI_MPICC=gcc-11


# COMPILER_VERSION := $(shell expr `gcc -dumpversion`)
ifeq ($(shell $(CC) -v 2>&1 | grep -c "clang version"), 1)
	XVERSION := $(shell cc -v 2>&1 | head -1 )
#	COMPILER_VERSION := $(shell $(CC) -v 2>&1 | head -1 )
else
#	COMPILER_VERSION := $(shell expr `gcc -dumpversion`)
	XVERSION:= gcc $(shell expr `gcc -dumpversion`)
endif


# This is mostly to address GitHub issues #100
# 221210 - Expression to find MPI Compiler modified so works with mpich as well as open mpi See issue #977
ifeq (mpicc, $(CC))
	# you're using mpicc, so we'll compile with the DMPI_ON flag
	MPI_FLAG = -DMPI_ON

	# check what underlying compiler mpi is using, and the version
	# we'll use this to print to the user and also to do some checks
	MPI_COMPILER := $(shell mpicc -show)
	COMPILER_VERSION := $(shell expr `'$(CC)' -dumpversion`)
	COMPILER_PRINT_STRING = Compiling with mpicc wrapper, for compiler $(MPI_COMPILER) $(COMPILER_VERSION)

	# if it's gcc we want to check if the version is 4.8 or later
	# if it is then we'll disable aggressive loop optimizations, see #100
	ifeq (gcc, $(MPI_COMPILER))
		GCCVERSIONGTEQ4 := $(shell expr `gcc -dumpversion | cut -f1-2 -d.` \>= 4.8)
		ifeq ("$(GCCVERSIONGTEQ4)", "1")
			EXTRA_FLAGS += -fno-aggressive-loop-optimizations 	# add the flag to EXTRA_FLAGS
			COMPILER_PRINT_STRING += with -fno-aggressive-loop-optimizations
		endif
	endif
else ifeq (gcc, $(CC))
	# no mpicc,
	MPI_FLAG =

	# check the version
	# we'll use this to print to the user and also to do some checks
	# COMPILER_VERSION := $(shell expr `gcc -dumpversion`)

	# if it's gcc we want to check if the version is 4.8 or later
	# if it is then we'll disable aggressive loop optimizations, see #100
	GCCVERSIONGTEQ4 := $(shell expr `gcc -dumpversion | cut -f1-2 -d.` \>= 4.8)
	COMPILER_PRINT_STRING = Compiling with $(CC) $(COMPILER_VERSION)
	ifeq ("$(GCCVERSIONGTEQ4)", "1")
		EXTRA_FLAGS += -fno-aggressive-loop-optimizations
		COMPILER_PRINT_STRING += with -fno-aggressive-loop-optimizations
	endif
else	# you must be using clang or icc
	MPI_FLAG =

	# check the version we'll use this to print to the user
	# COMPILER_VERSION = $(shell expr `$(CC) -dumpversion`)
	COMPILER_PRINT_STRING = Compiling with $(CC) $(COMPILER_VERSION)
endif

# If NVCC has been set, then we can define CUDA_FLAG = -DCUDA_ON, so the CUDA
# options in the C source are compiled
ifneq ($(NVCC), )
	CUDA_FLAG = -DCUDA_ON
	CUDA_SOURCE := $(wildcard *.cu)
	CUDA_LIBS = -lcudart -lcusolver
else
	CUDA_FLAG =
	CUDA_SOURCE =
	CUDA_LIBS =
endif

# These variables are used to compile CUDA code, and don't have to be defined
# conditionally
NVCC_FLAGS = -O3 -Werror all-warnings
CUDA_OBJECTS := $(patsubst %.cu,%.o,$(CUDA_SOURCE))

# this command finds out how many files with uncommitted changes there are
GIT_DIFF_STATUS := $(shell expr `git status --porcelain 2>/dev/null| grep "^ M" | wc -l`)
GIT_COMMIT_HASH := $(shell expr `git rev-parse HEAD`)

INDENT = yes
INCLUDE = ../include
LIB = ../lib
BIN = ../bin

ifeq (D,$(firstword $(MAKECMDGOALS)))
# use pg when you want to use gprof the profiler
# to use profiler make with arguments "make D sirocco"
# this can be altered to whatever is best
	CFLAGS = -std=gnu99 -g -pg -Wl,-Ttext-segment=0x68000000 -Wall -Werror $(EXTRA_FLAGS) -I$(INCLUDE) $(MPI_FLAG) ${CUDA_FLAG}
	FFLAGS = -g -pg
	PRINT_VAR = DEBUGGING, -g -pg -Wl,-Ttext-segment=0x68000000 -Wall flags
	XDEBUG = True
# Make the assumption that when using Clang the user is on MacOS, which doesn't
# have (easy?) access to the GNU profiler or CUDA
	ifeq ($(shell $(CC) -v 2>&1 | grep -c "clang version"), 1)
		CFLAGS = -std=gnu99 -g -Wall $(EXTRA_FLAGS) -I$(INCLUDE) $(MPI_FLAG)
		FFLAGS = -g
		PRINT_VAR = DEBUGGING, -g -Wall flags
	endif
else
# Use this for large runs
	CFLAGS = -std=gnu99 -O3 -Wall $(EXTRA_FLAGS) -I$(INCLUDE) $(MPI_FLAG) ${CUDA_FLAG}
	FFLAGS =
	PRINT_VAR = LARGE RUNS, -03 -Wall flags
endif



# next line for debugging when concerned about memory problems and duma installed in sirocco directory
# LDFLAGS= -L$(LIB)  -lm -lkpar  -lgslcblas ../duma_2_5_3/libduma.a -lpthread
# next line if you want to use kpar as a library, rather than as source below
# LDFLAGS= -L$(LIB)  -lm -lkpar -lcfitsio -lgsl -lgslcblas
LDFLAGS+= -L$(LIB) -lm -lgsl -lgslcblas $(CUDA_LIBS)


#Note that version should be a single string without spaces.


VERSION = 0.1



startup:
	@echo $(COMPILER_PRINT_STRING)			# prints out compiler information
	@echo 'Compiler is really  $(XVERSION)'
	@echo 'YOU ARE COMPILING FOR' $(PRINT_VAR)	# tells user if compiling for optimized or debug
	@echo 'CFLAGS='$(CFLAGS)
	@echo 'LDFLAGS='$(LDFLAGS)
	@echo 'MPI_FLAG='$(MPI_FLAG)
	@echo 'CUDA_FLAG='$(CUDA_FLAG)
	echo "#define VERSION " \"$(VERSION)\" > version.h
	echo "#define GIT_COMMIT_HASH" \"$(GIT_COMMIT_HASH)\" >> version.h
	echo "#define GIT_DIFF_STATUS" $(GIT_DIFF_STATUS)\ >> version.h

version:
	echo "#define VERSION " \"$(VERSION)\" > version.h
	echo "#define GIT_COMMIT_HASH" \"$(GIT_COMMIT_HASH)\" >> version.h
	echo "#define GIT_DIFF_STATUS" $(GIT_DIFF_STATUS)\ >> version.h

# To ensure a uniform style within .c programs, we enforce a indent style on .c code whe
# modifications are made. For this portion of the installation to work properly, one must have
# the python scripting installed, as well as gnuindent.  On linux systems, this is generally
# not an inssue, but this is not the case for MAC OSX.  The run_indent.py routine below tries
# to sort out whether one has gindent installed; if not it issues a woraning but does not
# reformat any of the changed files.

indent:
	# Fix indenting, if necessary.
	@if [ $(INDENT) = yes ] ; then  \
	../py_progs/run_indent.py -changed ; \
	fi \


zoo:
	# Fix indenting, if necessary.
	@if [ $(INDENT) = yes ] ; then  ../py_progs/run_indent.py -changed ; fi



foo: foo.o signal.o time.o
	$(CC) ${cfllags} foo.o signal.o time.o ${LDFLAGS}  -o foo

# For reasons that are unclear to me. get_models.c cannot be included in the sources.
# Problems occur due to the prototypes that are generated. Same for kpar_source it seems.
sirocco_source = agn.c anisowind.c atomic_extern_init.c atomicdata.c atomicdata_init.c  \
	atomicdata_sub.c bands.c bb.c bilinear.c brem.c cdf.c charge_exchange.c communicate_macro.c  \
	communicate_plasma.c communicate_spectra.c communicate_wind.c compton.c continuum.c cooling.c corona.c  \
	cv.c cylind_var.c cylindrical.c define_wind.c density.c diag.c dielectronic.c direct_ion.c  \
	disk.c disk_init.c disk_photon_gen.c emission.c estimators_macro.c estimators_simple.c  \
	extract.c frame.c  gradv.c gridwind.c homologous.c hydro_import.c import.c  \
	import_calloc.c import_cylindrical.c import_rtheta.c import_spherical.c ionization.c  \
	janitor.c knigge.c levels.c lines.c macro_accelerate.c macro_gen_f.c macro_gov.c  \
	matom.c matom_diag.c matrix_cpu.c matrix_ion.c models_extern_init.c para_update.c  \
	parse.c partition.c paths.c phot_util.c photon2d.c photon_gen.c photon_gen_matom.c  \
	pi_rates.c sirocco_extern_init.c radiation.c random.c rdpar.c rdpar_init.c recipes.c  \
	recomb.c resonate.c reverb.c roche.c rtheta.c run.c saha.c setup.c setup_disk.c setup_domains.c  \
	setup_files.c setup_line_transfer.c setup_reverb.c setup_star_bh.c shell_wind.c signal.c  \
	spectra.c spectral_estimators.c spherical.c stellar_wind.c sv.c synonyms.c time.c  \
	trans_phot.c vvector.c walls.c wind.c wind2d.c wind_sum.c wind_updates2d.c wind_util.c  \
	windsave.c windsave2table_sub.c xlog.c xtest.c zeta.c

# these are the objects required for compilation of sirocco. We are using pattern
# substitution so we don't have to maintain two identical lists but with .o instead of
# .c. Previously these two lists were slightly different. For some reason, get_models.c
# cannot be included in the sources, but we can define it OK in the objects.
sirocco_objects := $(sirocco_source:.c=.o) get_models.o

# kpar_source is now declared separately from sirocco_source so that the file log.h
# can be made using cproto
kpar_source = rdpar.c xlog.c synonyms.c
additional_swind_source = swind_sub.c swind_ion.c swind_write.c swind_macro.c swind.c windsave2table.c windsave2table_sub.c

prototypes:
	cp templates.h templates.h.old
	cproto -I$(INCLUDE) $(sirocco_source) ${additional_swind_source} > foo.h
	cp foo.h templates.h
	rm foo.h
	cproto -I$(INCLUDE) $(kpar_source) > log.h
	cproto  -I$(INCLUDE) atomicdata.c > atomic_proto.h
	cproto  -I$(INCLUDE) atomicdata_sub.c >> atomic_proto.h
	cproto  -I$(INCLUDE) atomicdata_init.c >> atomic_proto.h
	cproto  -I$(INCLUDE) recipes.c random.c cdf.c vvector.c > math_proto.h

# Recipe to create CUDA object code. If NVCC is blank, then nothing happens
# in this recipe
$(CUDA_OBJECTS): $(CUDA_SOURCE)
ifneq ($(CUDA_FLAG),)
	$(NVCC) $(NVCC_FLAGS) -DCUDA_ON -I$(INCLUDE) -c $< -o $@
endif

sirocco: startup sirocco.o $(sirocco_objects) $(CUDA_OBJECTS)
	$(CC) $(CFLAGS) sirocco.o $(sirocco_objects) $(CUDA_OBJECTS) $(kpar_objects) $(LDFLAGS) -o sirocco
	cp $@ $(BIN)/sirocco
        ifdef XDEBUG
	  mv $@ $(BIN)/sirocco-$(VERSION)_debug
        else
	  mv $@ $(BIN)/sirocco-$(VERSION)
        endif
	# Fix indenting, if necessary.
	@if [ $(INDENT) = yes ] ; then  ../py_progs/run_indent.py -changed ; fi


#This line is just so you can use make D sirocco for debugging
D:
	@echo 'Debugging Mode'

# Define the objects for the helper programs. They, in general, they all need
# to link with all of Python due to how interlinked a lot of it is.

swind_objects = swind.o swind_ion.o swind_macro.o swind_sub.o swind_write.o $(sirocco_objects)
table_objects = windsave2table.o $(sirocco_objects) # windsave2table_sub.o is already part of $(sirocco_objects)
rhf_objects = rad_hydro_files.o $(sirocco_objects)
modify_wind_objects = modify_wind.o $(sirocco_objects)
inspect_wind_objects = inspect_wind.o $(sirocco_objects)
sirocco_optd_obj = sirocco_optd.o sirocco_optd_output.o sirocco_optd_trans.o sirocco_optd_util.o sirocco_optd_extern_init.o $(sirocco_objects)
windsave2fits_objects = windsave2fits.o $(sirocco_objects)

run_indent:
	../py_progs/run_indent.py -all_no_headers

swind: startup $(swind_objects) $(CUDA_OBJECTS)
	$(CC) $(CFLAGS) $(swind_objects) $(CUDA_OBJECTS) $(LDFLAGS) -o swind
	cp $@ $(BIN)
	mv $@ $(BIN)/swind$(VERSION)
	@if [ $(INDENT) = yes ] ; then  ../py_progs/run_indent.py -changed ; fi

windsave2table: startup $(table_objects) $(CUDA_OBJECTS)
	$(CC) $(CFLAGS) $(table_objects) $(CUDA_OBJECTS) $(LDFLAGS) -o windsave2table
	cp $@ $(BIN)
	mv $@ $(BIN)/windsave2table$(VERSION)
	@if [ $(INDENT) = yes ] ; then  ../py_progs/run_indent.py -changed ; fi

windsave2fits: startup $(windsave2fits_objects) $(CUDA_OBJECTS)
	$(CC) $(CFLAGS) $(windsave2fits_objects) $(CUDA_OBJECTS) $(LDFLAGS) -lcfitsio -o windsave2fits 
	cp $@ $(BIN)
	mv $@ $(BIN)/windsave2fits$(VERSION)
	@if [ $(INDENT) = yes ] ; then  ../py_progs/run_indent.py -changed ; fi

rad_hydro_files: startup $(rhf_objects) $(CUDA_OBJECTS)
	$(CC) $(CFLAGS) $(rhf_objects) $(CUDA_OBJECTS) $(LDFLAGS) -o rad_hydro_files
	cp $@ $(BIN)
	mv $@ $(BIN)/rad_hydro_files$(VERSION)
	@if [ $(INDENT) = yes ] ; then  ../py_progs/run_indent.py -changed ; fi

modify_wind: startup $(modify_wind_objects) $(CUDA_OBJECTS)
	$(CC) $(CFLAGS) $(modify_wind_objects) $(CUDA_OBJECTS) $(LDFLAGS) -o modify_wind
	cp $@ $(BIN)
	mv $@ $(BIN)/modify_wind$(VERSION)
	@if [ $(INDENT) = yes ] ; then  ../py_progs/run_indent.py -changed ; fi

inspect_wind: startup $(inspect_wind_objects) $(CUDA_OBJECTS)
	$(CC) $(CFLAGS) $(inspect_wind_objects) $(CUDA_OBJECTS) $(LDFLAGS)  -o inspect_wind
	cp $@ $(BIN)
	mv $@ $(BIN)/inspect_wind$(VERSION)
	@if [ $(INDENT) = yes ] ; then  ../py_progs/run_indent.py -changed ; fi

sirocco_optd: startup $(sirocco_optd_obj) $(CUDA_OBJECTS)
	$(CC) $(CFLAGS) $(sirocco_optd_obj) $(CUDA_OBJECTS) $(LDFLAGS) -o $@
	cp $@ $(BIN)
	mv $@ $(BIN)/sirocco_optd$(VERSION)
	@if [ $(INDENT) = yes ] ; then  ../py_progs/run_indent.py -changed ; fi

# unit_tests allows one to test low level routines, but one must be aware that one needs to define
# everything in the main routine before one can probe a particular routine.
unit_test: startup unit_test.o $(sirocco_objects)
	$(CC) ${CFLAGS} unit_test.o $(sirocco_objects) $(kpar_objects) $(LDFLAGS) -o unit_test
	cp $@ $(BIN)/unit_test
	# Fix indenting, if necessary.
	@if [ $(INDENT) = yes ] ; then  ../py_progs/run_indent.py -changed ; fi

# The next recipe is used to run the unit tests, using the same CC and NVCC as
# the Python source. It doesn't require anything to be compiled, actually, as
# the Makefile in tests/ will compile what is requried
check:
	cd tests; make check NVCC=$(NVCC) CC=$(CC) MATOM_VER=$(MATOM_VER)

# The next line runs recompiles all of the routines after first cleaning the directory
# all: clean run_indent sirocco windsave2table swind
all: clean sirocco windsave2table swind indent rad_hydro_files modify_wind inspect_wind sirocco_optd

FILE = atomicdata.o atomicdata_init.o atomicdata_sub.o atomic.o

libatomic.a:  atomicdata.o atomicdata_init.o atomic.o
	aor ru libatomic.a atomicdata.o atomic.o
	ranlib libatomic.a
	mv libatomic.a $(LIB)
	cp atomic.h  $(INCLUDE)

clean :
	rm -f *.o  *~
	cd tests; make clean
