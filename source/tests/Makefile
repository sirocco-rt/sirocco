C = gcc
NVCC =
MATOM_VER = 1
TARGET = sirocco_unit_test

I_FLAGS = -I$(SIROCCO)/include
L_FLAGS = -L$(SIROCCO)/lib
C_FLAGS = -O3 -Wall -Wno-deprecated-non-prototype -DMATOM_VER=$(MATOM_VER)
CU_FLAGS = -O3

INCLUDES =
LIBS = -lstdc++ -lcunit -lgsl -lgslcblas -lm

# This is less than ideal, but it's far too complicated and time consuming to
# figure out what tests needs which source files. So as we've done in the past,
# we'll just include ALL of Python's source
SIROCCO_SOURCES = \
	agn.c anisowind.c atomic_extern_init.c atomicdata.c atomicdata_init.c  \
	atomicdata_sub.c bands.c bb.c bilinear.c brem.c cdf.c charge_exchange.c  \
	compton.c continuum.c cooling.c corona.c cv.c cylind_var.c cylindrical.c  \
	density.c diag.c dielectronic.c direct_ion.c disk.c disk_init.c  \
	disk_photon_gen.c emission.c estimators_macro.c estimators_simple.c  \
	extract.c frame.c gradv.c gridwind.c homologous.c hydro_import.c import.c  \
	import_calloc.c import_cylindrical.c import_rtheta.c import_spherical.c  \
	ionization.c knigge.c levels.c lines.c macro_accelerate.c macro_gen_f.c  \
	macro_gov.c matom.c matom_diag.c matrix_cpu.c matrix_ion.c  \
	models_extern_init.c para_update.c parse.c partition.c paths.c  \
	phot_util.c photon2d.c photon_gen.c photon_gen_matom.c pi_rates.c  \
	sirocco_extern_init.c radiation.c random.c rdpar.c recipes.c recomb.c  \
	resonate.c reverb.c roche.c rtheta.c run.c saha.c setup.c setup_disk.c  \
	setup_domains.c setup_files.c setup_line_transfer.c setup_reverb.c  \
	setup_star_bh.c shell_wind.c signal.c spectra.c spectral_estimators.c  \
	spherical.c stellar_wind.c sv.c synonyms.c time.c trans_phot.c  \
	vvector.c walls.c wind.c wind2d.c wind_sum.c wind_updates2d.c wind_util.c  \
	windsave.c windsave2table_sub.c xlog.c xtest.c zeta.c get_models.c \
	rdpar_init.c define_wind.c communicate_wind.c communicate_plasma.c \
	communicate_macro.c communicate_spectra.c

# This is sources to be compiled with NVCC
CUDA_SOURCES = \
	matrix_gpu.cu

# Source code for unit tests
TEST_SOURCES = \
	unit_test_main.c \
	unit_test_model.c \
	tests/test_matrix.c \
	tests/test_compton.c \
	tests/test_define_wind.c \
	tests/test_run_mode.c \
	tests/test_translate.c

# Using absolute paths
SIROCCO_SOURCES := $(patsubst %,$(SIROCCO)/source/%, $(SIROCCO_SOURCES))
CUDA_SOURCES := $(patsubst %,$(SIROCCO)/source/%, $(CUDA_SOURCES))
TEST_SOURCES := $(patsubst %,$(SIROCCO)/source/tests/%, $(TEST_SOURCES))

# This is all sources, to be compiled with CC
C_SOURCES := $(TEST_SOURCES) $(SIROCCO_SOURCES)

# Object files, again keeping C and CUDA sources separate
C_OBJS := $(C_SOURCES:.c=.o)

ifeq ($(CC), mpicc)
	C_FLAGS += -DMPI_ON
endif

ifneq ($(NVCC),)
	CUDA_OBJS := $(CUDA_SOURCES:.cu=.o)
	LIBS += -lcudart -lcusolver
	C_FLAGS += -DCUDA_ON
	CU_FLAGS += -DCUDA_ON
else
	CUDA_OBJS =
endif

all: startup clean $(TARGET)

startup:
	@echo "CC        =" $(CC)
	@echo "NVCC      =" $(NVCC)
	@echo "I_FLAGS   =" $(I_FLAGS)
	@echo "L_FLAGS   =" $(L_FLAGS)
	@echo "C_FLAGS   =" $(C_FLAGS)
	@echo "CU_FLAGS  =" $(CU_FLAGS)
	@echo "INCLUDES  =" $(INCLUDES)
	@echo "LIBS      =" $(LIBS)

$(CUDA_OBJS): %.o : %.cu
ifneq ($(NVCC),)
	$(NVCC) $(CU_FLAGS) $(I_FLAGS) $(INCLUDES) -c $< -o $@
endif

$(C_OBJS): %.o : %.c
	$(CC) $(C_FLAGS) $(I_FLAGS) $(INCLUDES) -c $< -o $@

$(TARGET): $(C_OBJS) $(CUDA_OBJS)
	$(CC) $(C_OBJS) $(CUDA_OBJS) $(L_FLAGS) $(LIBS) -o $@
	mv $@ $(SIROCCO)/bin

check: startup $(TARGET)
	$(SIROCCO)/bin/$(TARGET)

clean:
	rm -f *.o
	rm -f tests/*.o
	rm -f ../*.o

.PHONY: clean
