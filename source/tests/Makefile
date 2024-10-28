include $(SIROCCO)/source/Makefile 
# This is less than ideal, but it's far too complicated and time consuming to
# figure out what tests needs which source files. So as we've done in the past,
# we'll just include ALL of Python's source. 
# This is done by including the sirocco_source from the main Makefile 
# note that for some reason, get_models.c cannot be included in the sources
# in that Makefile, so we have to add it in
SIROCCO_SOURCES = $(sirocco_source) get_models.c


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
