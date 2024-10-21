Matrix Acceleration using CUDA
###############################

SIROCCO can use CUDA/GPUs to accelerate solving linear systems using the cuSOLVER library, which is part of the NVIDIA
CUDA Toolkit.

*This pilot study into using GPUs in SIROCCO was conducted as an HPC RSE project at the University of Southampton.*

When should you use CUDA?
=========================

Given the pilot study nature of this work, the current matrix acceleration implementation (September 2023) is simple.
In most small to mid-sized models, using GPU acceleration will not improve model performance. As of writing,
CUDA is only used to accelerate matrix calculations specifically, so there are no performance improvements other than in
models which are bogged down by matrix calculations; e.g. such as when calculating the ionization state for a large
number of ions, or models with lots of macro atom levels/emissitivies. Even so, there may only be modest improvements in
model performance if the majority of time is spent transporting and scattering photons.

It is therefore only preferable to use CUDA when matrices are large enough to warrant GPU acceleration. The size of
where this is important is tricky to know, as it is hardware dependent - both on your CPU and GPU. If you are using an
old CPU, then you are likely to see improvements from matrices as small as circa 200 x 200. Otherwise, you may need to
reach up to matrix sizes of 500 x 500 (or larger!) before there is any tangible benefit.

This situation will be improved with smarter memory management and further parallelisation. Most time is spent copying
data back and forth between the GPU and CPU. As an example, consider the matrix ionization state calculation. Currently
only the actual step to stop the linear system (to calculate the ionization state) has been ported to the GPU. This
means each iteration toward a converged ionization state requires memory to be copied to and from the GPU, which slows
things down quite a bit. If you could instead port the entire iterative procedure to the GPU (which is not that easy),
there is no longer the need to make expensive memory copies each iteration which will significantly speed up the
algorithm.

Requirements
============

To use the CUDA matrix acceleration in SIROCCO, your machine needs to have the following installed,

- A CUDA-capable NVIDIA GPU
- NVIDIA CUDA Toolkit
- NVIDIA GPU drivers
- A supported operating system (Windows or Linux) with a gcc compiler and toolchain

NVIDIA provides a list of CUDA-enabled GPUs `here <https://developer.nvidia.com/cuda-gpus>`__. Whilst the GeForce series
of NVIDIA GPUs are more affordable and generally *good enough*, from a purely raw computation standpoint NVIDIA's
workstation and data center GPUs are more suited due differences (and additional) in hardware not included in
the GeForce line of GPUs.

Installing the CUDA toolkit
---------------------------

The NVIDIA CUDA Toolkit is installed either through an installer downloaded from NVIDIA or can be installed via a
package manager on Linux systems. It should be noted that the CUDA Toolkit *does not* come with NVIDIA drivers and need
to be installed separately. The NVIDIA CUDA Toolkit is available at `https://developer.nvidia.com/cuda-downloads
<https://developer.nvidia.com/cuda-downloads>`_ and NVIDIA drivers at `https://www.nvidia.co.uk/Download/index.aspx
<https://www.nvidia.co.uk/Download/index.aspx?lang=en-uk>`_.

On Ubuntu 22.04, the toolkit and NVIDIA's proprietary drivers are available through :code:`apt`,

.. code:: bash

    sudo apt install nvidia-cuda-toolkit nvidia-driver-535

How to Enable and Run CUDA
==========================

Compilation
-----------

CUDA is an additional acceleration method and is therefore not enabled by default. To enable CUDA, SIROCCO has to be
compiled with the additional :code:`-DCUDA_ON` flag and linked with the appropriate libraries using the NVIDIDA CUDA
compiler (nvcc). There are several ways to enable the CUDA components of SIROCCO. The most simple is to run the configure
script in the root SIROCCO directory with the arguments :code:`--with-cuda`,

.. code:: bash

    [$SIROCCO] $ ./configure --with-cuda

    Configuring Makefile for SIROCCO radiative transfer code
    Checking for mpicc...yes
    Checking for gcc...yes
    Checking for nvcc...yes
    Preparing Makefile...Done.

If the NVIDIA CUDA Toolkit is found, you will see the output informing that the CUDA compiler :code:`nvcc` was found.

What essentially happens when you run :code:`code` is that a value for the variable :code:`NVCC` is set in the Makefile
in SIROCCO's source directory. If you re-run :code:`configure` without :code:`--with-cuda`, then :code:`NVCC` will be
unset and CUDA will not be used. CUDA can be disabled or enabled *"on the fly"* by modifying this variable without
running the configure script and by modifying the Makefile or passing the value of the variable when calling the
Makefile, e.g.,

.. code:: bash

    [$SIROCCO/source] $ make clean
    [$SIROCCO/source] $ make sirocco NVCC=nvcc

:code:`make clean` has to be run whenever CUDA is switched been enabled or disabled, due to code conditionally compiling
depending on if CUDA is enabled or not.

Running
-------

To run SIROCCO with CUDA, you run it in the exact way even parallel models running with MPI. On a HPC system the
appropriate GPU resources will need to be requested in a job submission script. For example, on Iridis at the University
of Southampton, a functional job submission script may look like this,

.. code:: bash

    #!/bin/bash

    #SBATCH --nodes=1
    #SBATCH --ntasks-per-node=40
    #SBATCH --gpus-per-node=1
    #SBATCH --time=06:00:00
    #SBATCH --partition=gpu

    module load openmpi/4.1.5/gcc

    mpirun -n $SLURM_NTASKS py model.pf

If CUDA is enabled and no GPU resources are found, SIROCCO will exit early in the program with an appropriate error
message. Note that a CUDA-aware MPI implementation is not required, as no data is communicated between GPUs.

Implementation
==============

In this part of the documentation, we will cover the implementation details of cuSolver in SIROCCO. cuSolver is a matrix
library within the NVIDIA CUDA ecosystem, designed to accelerate both dense and sparse linear algebra problems,
including matrix factorisation, linear system solving and matrix inversion. To use cuSolver, very little GPU specific
code needs to be written, other than code to allocate memory on the GPU. There are therefore a number of similarities
between writing functions which use the cuSolver (and other CUDA mathematical libraries) and GSL libraries.

The CUDA parallel model
-----------------------

The main difference between CPU and GPU parallel programming is the number of (dumb) cores in a GPU. Whereas on a CPU
where we divide work on a matrix into smaller chunks, on a GPU it is realistic to have each core of the GPU operate on a
single element of the matrix whereas a CPU will likely have multiple elements. CUDA is a type of shared memory parallel
programming, and at its core are kernels, which are specialised functions designed for massive parallelism. These
kernels are executed by each thread (organized in blocks and grids), where thousands are launched and execute the code
concurrently allowing for massive parallelism.

As an example, consider matrix multiplication. If the calculation is parallelised, each CPU core will likely need to
calculate the matrix product for multiple elements of the matrix. On a GPU, each thread that is launched will calculate
the product for only a single element. If there are enough GPU cores available, then the calculation can be done in
effectively a single step which all threads calculating the product for each element at once.

A more detailed and thorough explanation of the CUDA programming model can be found in the `CUDA documentation
<https://docs.nvidia.com/cuda/cuda-c-programming-guide/index.html#>`__.

Basics
------

SIROCCO uses the dense matrix functions in cuSolver, which are accessed through the :code:`cusolverDn.h` header file. To
use cuSolver, it must first be initialised. To do so, we use :code:`cusolverDnCreate` to create a
:code:`cuSolverDnHandle_t` variable which is used by cuSolver internally for resource and context management.

cuSolver is based on the Fortran library `LAPACK <https://www.netlib.org/lapack/>`_ and as such expects arrays to be
ordered in column-major order like in Fortran. In C, arrays are typically ordered in row-major order and so arrays must
be transposed into column-major ordering before being passed to cuSolver (an explanation of the differences between row
and column major ordering can be found `here <https://en.wikipedia.org/wiki/Row-_and_column-major_order>`__). Matrices
can be transposed either whilst still on the CPU, or on the GPU by using a CUDA kernel as shown in the example below,

.. code:: cpp
    :caption: A CUDA kernel to transpose a matrix from row to column major

    __global__ void  /* __global__ is used by kernels, all of which return void */
    transpose_row_to_column_major(double *row_major, double *column_major, int matrix_size)
    {
        /* Determine the x and y coordinate for the thread -- these coords could be
           outside the matrix if enough threads are spawned */
        const int idx = blockIdx.x * blockDim.x + threadIdx.x;
        const int idy = blockIdx.y * blockDim.y + threadIdx.y;

        /* Only transpose for threads inside the matrix */
        if (idx < matrix_size && idy < matrix_size) {
            column_major[idx * matrix_size + idy] = row_major[idy * matrix_size + idx];
        }
    }

The syntax of the above is covered in detail in the `CUDA documentation
<https://docs.nvidia.com/cuda/cuda-c-programming-guide/index.html#kernels>`__. The purpose of the kernel is take in a row
major array and to transpose it to column major.

Every cuSolver (and CUDA) function returns an error status. To make code more readable, a macro is usually defined which
checks the error status and raises an error message if the function does not execute successfully. This type of macro is
used extensively throughout the implementation.

.. code:: c
    :caption: A useful macro for error checking cuSolver returns

    #define CUSOLVER_CHECK(status)                                                                                     \
        do {                                                                                                           \
            cusolverStatus_t err = status;                                                                             \
            if (err != CUSOLVER_STATUS_SUCCESS) {                                                                      \
                Error("cuSolver Error (%d): %s (%s:%d)\n", err, cusolver_get_error_string(err), __FILE__, __LINE__);   \
                return err;                                                                                            \
            }                                                                                                          \
        } while (0)

    /* Here is an example of using the macro to create a handle */
    CUSOLVER_CHECK(cusolverDnCreate(&handle));

Structure
---------

When writing CUDA C, it is convention to put the CUDA code into :code:`.cu` files and the CPU code in :code:`.c` files.
Even when using a library like cuSolver, it is still convention to place that code into :code:`.cu` files as we still
need to access some CUDA library functions, such as :code:`cudaMalloc` or :code:`cudaMemCpy`.

The CUDA code associated with matrix parallelisation has been written in the file :code:`$SIROCCO/source/matrix_gpu.cu`
with the header file :code:`$SIROCCO/source/matrix_gpu.h` which includes the function prototypes for the GPU matrix code.
The GSL matrix code is kept in :code:`$SIROCCO/source/matrix_cpu.c` with function prototypes in
:code:`$SIROCCO/source/templates.h`.

To be able to switch between the CUDA and GSL matrix implementations with the minimal amount of code changes, a
:code:`solve_matrix` wrapper function has been created. Either GSL or cuSolver is called within this wrapper, depending
on if SIROCCO was compiled with the flag :code:`-DCUDA_ON` as discussed earlier. This wrapper takes on the same name as
the original GSL implementation, meaning no code changes have occurred in that regard.

.. code:: c
    :caption: The wrapper function which calls the appropriate matrix solver

    #include "matrix_gpu.h"  /* The function prototype for gpu_solve_matrix is in here */

    int
    solve_matrix(double *a_matrix, double *b_vector, int matrix_size, double *x_vector)
    {
        int error;

    #ifdef CUDA_ON
        error = gpu_solve_matrix(...);  /* CUDA implementation */
    #else
        error = cpu_solve_matrix(...);  /* GSL implementation */
    #endif

        return error;
    }

The following code exert is an example of using the wrapper function to solve a linear system.

.. code:: c
    :caption: The API to solve a linear system hasn't changed

    #include "python.h"

    double *populations = malloc(nions * sizeof(*populations));
    double *ion_density = malloc(nions * sizeof(*populations));
    double *rate_matrix = malloc(nions * nions * sizeof(*populations));

    populate_matrices(rate_matrix, ion_density);

    /* The wrapper function is named the same as the original GSL implementation
       and accepts the same arguments */
    int error = solve_matrix(
        rate_matrix, ion_density, nions, populations, xplasma->nplasma
    );

    /* One user difference is that error handling is more robust now, and there
       is a function to convert error codes into error messages */
    if (error != EXIT_SUCCESS) {
        Error(
            "Error whilst solving for ion populations: %d (%d)\n",
            get_matrix_error_string(error), error
        );
    }

Here is an example of using a similar wrapper function to calculate the inverse of a matrix.

.. code:: c
    :caption: The API has changed slightly for calculating the inverse, now that it has a wrapper function

    #include "python.h"

    double Q_matrix = malloc(matrix_size * matrix_size * sizeof(double));
    double Q_inverse = malloc(matrix_size * matrix_size * sizeof(double));

    populate_Q_matrix(Q_matrix);

    /* The API is only different in the sense that a wrapper function now
       exists for matrix inversion */
    int error = invert_matrix(
        Q_matrix, Q_inverse, matrix_size
    );

    if (error != EXIT_SUCCESS) {
        Error(
            "Error whilst solving for ion populations: %d (%d)\n",
            get_matrix_error_string(error), error
        );
    }

To write the cuSolver implementation is similar to the GSL implementation, in that memory/resource are allocated for
cuSolver and then the appropriate library functions are called. The code exert below shows an illustrated (and
simplified) example of the cuSolver implementation to solve a linear system.

.. code:: c
    :caption: An illustrative example of using cuSolver to solve a linear system using LU decomposition

    #include <cuSolverDn.h>

    extern "C" int  /* extern "C" has to be used to make it available to the C run time */
    gpu_solve_matrix(double *a_matrix, double *b_vector, int matrix_size, double *x_vector)
    {
        /* First of all, allocate memory on the GPU and copy data from the CPU to the
           GPU. This uses the CUDA standard library functions, such as cudaMemCpy and
           cudaMalloc. This is part of the code is what takes the most time. */
        allocate_memory_for_gpu();
        copy_data_to_gpu();

        /* cuSolver and cuBLAS are both ports of Fortran libraries, which expect arrays to
        be in column-major format and we therefore need to transpose our row-major arrays */
        transpose_row_to_column_major<<<grid_dim, block_dim>>>(
            d_matrix_row, d_matrix_col, matrix_size
        );

        /* Perform LU decomposition. Variables prefixed with d_ are kept in GPU memory where we
        allocated space for them in `allocate_memory_for_gpu` */
        CUSOLVER_CHECK(cusolverDnDgetrf(
            CUSOLVER_HANDLE, matrix_size, matrix_size, d_matrix_col, matrix_size,
            d_workspace, d_pivot, d_info
        ));

        /* Solve the linear system A x = b. The final solution is returned in the
           variable d_v_vector */
        CUSOLVER_CHECK(cusolverDnDgetrs(
            CUSOLVER_HANDLE, CUSOLVER_OP_N, matrix_size, matrix_size, d_matrix_col,
            matrix_size, d_pivot,
            d_b_vector, matrix_size, d_info
        ));

        /* We now have to copy d_b_vector back to the CPU, so we can use that value in
        the rest of SIROCCO */
        copy_data_to_cpu();

        return EXIT_SUCCESS;
    }

The naming conventions of cuSolver are discussed `here
<https://docs.nvidia.com/cuda/cusolver/index.html#naming-conventions>`__. In the case above, :code:`cuSolverDnDgetrf`
corresponds to: cusolverDn = *cuSolver Dense Matrix*, D = *double precision (double)* and getrf = *get right
hand factorisation*.

The most important thing to note, which may appear trivial, is the :code:`extern` keyword. Without this, when the
program is compiled the function :code:`gpu_solve_matrix` will not be available to the C runtime. By labelling the
function as :code:`extern "C"`, we make it clear that we want this function to be available to C source code. This only
needs to be done at the function definition, and not the function prototype in, e.g., a header file.

Compiling and Linking
---------------------


CUDA code is compiled using the NVIDIA CUDA Compiler :code:`nvcc`. To combine both CPU and GPU code, the source must be
compiled with the respective compilers (e.g. :code:`gcc`/:code:`mpicc` for C and :code:`nvcc` for CUDA) to object code
(:code:`.o` files) and which are linked together using the C compiler with the appropriate library flags. In addition to
needing to link the cuSolver library (:code:`-lcusolver`) we also need to link the CUDA runtime library
(:code:`-lcudart`) when linking with the C compiler, which makes the standard CUDA library functions available to the C
compiler and runtime.

The steps for compiling and link GPU and CPU code are outlined below in pseudo-Makefile code.

.. code:: bash
    :caption: A brief overview on how to compile and link C and CUDA code

    # Define compilers for C and CUDA. When creating a CUDA/MPI application, we can
    # just as easily use mpicc for our C compiler. It makes no difference.
    CC = mpicc
    NVCC = nvcc

    # Define C and CUDA libraries. We still include GSL as other GSL numerical routines are
    # used in SIROCCO
    C_LIBS = -lgsl -lgslcblas -lm
    CUDA_LIBS = -lcudart -lcusolver

    # Define flags for C and CUDA compilers. -DCUDA_ON is used to conditionally compile
    # to use the CUDA wrappers and other things related to the CUDA build
    C_FLAGS = -O3 -DCUDA_ON -DMPI_ON -I../includes -L../libs
    CUDA_FLAGS = -O3 -DCUDA_ON

    # Compile CUDA source to object code using the CUDA compiler
    $(NVCC) $(CUDA_FLAGS) $(CUDA_SOURCE) -c -o $(CUDA_OBJECTS)

    # Compile the C code using the C compiler
    $(CC) $(C_FLAGS) $(C_SOURCE) -c -o $(C_OBJECTS)

    # Link the CUDA and C object code and libraries together using the C compiler
    $(CC) $(CUDA_OBJECTS) $(C_OBJECTS) -o python $(CUDA_LIBS) $(C_LIBS)

These steps are effectively replicated in the Makefile :code:`$SIROCCO/source/Makefile`, where a deconstructed example is
shown below.

.. code:: Makefile
    :caption: The variables and recipes associated with CUDA are all conditional on NVCC being defined

    # If NVCC has been set in the Makefile, then we can define CUDA_FLAG = -DCUDA_ON,
    # and the CUDA sources, which, at the moment, uses a wildcard to find all .cu files
    ifneq ($(NVCC), )
    	CUDA_FLAG = -DCUDA_ON
    	CUDA_SOURCE = $(wildcard *.cu)
    else
    	CUDA_FLAG =
    	CUDA_SOURCE =
    endif

    # Then the recipe to create CUDA object code looks like this. If NVCC is blank,
    # nothing happens in the recipe
    $(CUDA_OBJECTS): $(CUDA_SOURCE)
    ifneq ($(CUDA_FLAG),)
        $(NVCC) $(NVCC_FLAGS) -DCUDA_ON -I$(INCLUDE) -c $< -o $@
    endif

    # So to compile SIROCCO, we have something which looks vaguely like this. Note that
    # we use the CUDA_OBJECTS recipe as a requirement for the python recipe. This CUSOLVER_STATUS_SUCCESS
    # the CUDA source to be compiled to object code *if* NVCC is defined
    python: startup python.o $(python_objects) $(CUDA_OBJECTS)
        $(CC) $(CFLAGS) python.o $(python_objects) $(CUDA_OBJECTS) $(kpar_objects) $(LDFLAGS) -o python

