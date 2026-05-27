# Python CMakeLists
#
# The main purpose of this CMake file is to be able to build some of the
# binaries part of the Python project. It can't handle building from scratch, so
# you will probably need to build things using the Makefile first of all.
#
# So really, this is only helpful if you need to define things for a CMake
# project in VSCode or CLion, for example. It's good enough to define programs,
# and build for debugging. But isn't complex enough to fully emulate the
# behaviour of the Makefile in $PYTHON/source.

cmake_minimum_required(VERSION 3.5)
project(Python LANGUAGES C)

# The C standard is being set to 90. The compiler is being set to mpicc

set(CMAKE_C_STANDARD 99)
set(CMAKE_C_COMPILER mpicc)
set(CMAKE_CXX_COMPILER mpixx)
include_directories(include)
link_directories(lib)
add_definitions(-DMPI_ON)

# Set the Python version, which is going to be appended to source/version.h

set(PYTHON_VERSION 87g)
message(STATUS "Python version : ${PYTHON_VERSION}")

# Create the version.h file. During the configure process, the file version.h
# should be created and whenever the program is build again, the information in
# version.h should be updated.

set(GIT_HASH "")
set(GIT_NUM_CHANGES 0)
set(VERSION_FILE source/version.h)

if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/.git)
  find_package(Git)
  find_program(BASH bash HINTS /bin)

  if(GIT_FOUND)
    execute_process(
      COMMAND ${GIT_EXECUTABLE} rev-parse HEAD
      WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
      OUTPUT_VARIABLE GIT_HASH
      ERROR_QUIET OUTPUT_STRIP_TRAILING_WHITESPACE)
    execute_process(
      COMMAND ${BASH} "-c"
              "${GIT_EXECUTABLE} status --porcelain | grep '^ M' | wc -l"
      WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
      OUTPUT_VARIABLE GIT_NUM_CHANGES
      ERROR_QUIET OUTPUT_STRIP_TRAILING_WHITESPACE)
  else(GIT_FOUND)
    set(GIT_HASH 0)
    set(GIT_NUM_CHANGES 0)
  endif(GIT_FOUND)

  message(STATUS "Git hash       : ${GIT_HASH}")
  message(STATUS "Git changes    : ${GIT_NUM_CHANGES}")

  file(
    WRITE ${VERSION_FILE}
    "#define VERSION \"${PYTHON_VERSION}\"\n\
#define GIT_COMMIT_HASH \"${GIT_HASH}\"\n\
#define GIT_DIFF_STATUS \"${GIT_NUM_CHANGES}\"\n")
endif(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/.git)

# Define the source for the various different programs

set(CUDA_SOURCE source/matrix_gpu.cu)

set(PYTHON_SOURCE
    ${CUDA_SOURCE}
    source/agn.c
    source/atomicdata.c
    source/atomicdata_init.c
    source/atomicdata_sub.c
    source/anisowind.c
    source/bands.c
    source/bb.c
    source/bilinear.c
    source/brem.c
    source/cdf.c
    source/charge_exchange.c
    source/compton.c
    source/continuum.c
    source/cooling.c
    source/corona.c
    source/cv.c
    source/cylind_var.c
    source/cylindrical.c
    source/density.c
    source/diag.c
    source/dielectronic.c
    source/direct_ion.c
    source/disk.c
    source/disk_init.c
    source/emission.c
    source/estimators_macro.c
    source/estimators_simple.c
    source/extract.c
    source/frame.c
    source/get_models.c
    source/gradv.c
    source/gridwind.c
    source/homologous.c
    source/hydro_import.c
    source/import.c
    source/import_calloc.c
    source/import_cylindrical.c
    source/import_polar.c
    source/import_spherical.c
    source/ionization.c
    source/knigge.c
    source/levels.c
    source/lines.c
    source/macro_gov.c
    source/matom.c
    source/matom_diag.c
    source/matrix_ion.c
    source/para_update.c
    source/parse.c
    source/partition.c
    source/paths.c
    source/phot_util.c
    source/photon2d.c
    source/photon_gen.c
    source/pi_rates.c
    source/radiation.c
    source/random.c
    source/rdpar.c
    source/rdpar_init.c
    source/recipes.c
    source/recomb.c
    source/resonate.c
    source/reverb.c
    source/roche.c
    source/rtheta.c
    source/run.c
    source/saha.c
    source/setup.c
    source/setup_disk.c
    source/setup_domains.c
    source/setup_files.c
    source/setup_line_transfer.c
    source/setup_reverb.c
    source/setup_star_bh.c
    source/shell_wind.c
    source/signal.c
    source/spectra.c
    source/spectral_estimators.c
    source/spherical.c
    source/stellar_wind.c
    source/sv.c
    source/synonyms.c
    source/time.c
    source/trans_phot.c
    source/vvector.c
    source/walls.c
    source/wind.c
    source/wind2d.c
    source/wind_sum.c
    source/wind_updates2d.c
    source/wind_util.c
    source/windsave.c
    source/windsave2table_sub.c
    source/xlog.c
    source/xtest.c
    source/zeta.c
    source/macro_accelerate.c
    source/macro_gen_f.c
    source/atomic_extern_init.c
    source/python_extern_init.c
    source/models_extern_init.c
    source/photon_gen_matom.c
    source/disk_photon_gen.c
    source/matrix_cpu.c
    source/define_wind.c
    source/communicate_wind.c
    source/communicate_plasma.c
    source/communicate_macro.c
    source/communicate_spectra.c
    source/janitor.c
)

set(PY_WIND_SOURCE source/py_wind_ion.c source/py_wind_macro.c
                   source/py_wind_sub.c source/py_wind_write.c
        source/define_wind.c
        source/communicate_wind.c
        source/communicate_spectra.c)

set(WINDSAVE2TABLE_SOURCE source/windsave2table_sub.c
        source/define_wind.c
        source/communicate_wind.c)

set(RAD_HYDRO_SOURCE source/rad_hydro_files.c)

set(MODIFY_WIND_SOURCE source/modify_wind.c)

set(INSPECT_WIND_SOURCE source/inspect_wind.c)

set(OPTICAL_DEPTH_SOURCE
    source/py_optd.c source/py_optd_trans.c source/py_optd_util.c
    source/py_optd_output.c source/py_optd_extern_init.c)

set(TEST_SOURCE
    source/tests/unit_test_main.c
    source/tests/tests/test_matrix.c
    source/tests/tests/test_compton.c
    ${PYTHON_SOURCE}
        source/define_wind.c
        source/tests/tests/test_define_wind.c
        source/communicate_wind.c
        source/tests/unit_test_model.c
        source/tests/tests/test_run_mode.c
)

# Create the executables for each program, and link the required libraries

# Python
add_executable(python source/python.c ${PYTHON_SOURCE})
target_link_libraries(python gsl gslcblas m)

# Py_wind
add_executable(py_wind source/py_wind.c ${PYTHON_SOURCE} ${PY_WIND_SOURCE})
target_link_libraries(py_wind gsl gslcblas m)

# Windsave2table
add_executable(windsave2table source/windsave2table.c ${PYTHON_SOURCE}
                              ${WINDSAVE2TABLE_SOURCE})
target_link_libraries(windsave2table gsl gslcblas m)

# rad_hydro_files
add_executable(rad_hydro_files ${PYTHON_SOURCE} ${RAD_HYDRO_SOURCE})
target_link_libraries(rad_hydro_files gsl gslcblas m)

# modify_wind
add_executable(modify_wind ${PYTHON_SOURCE} ${MODIFY_WIND_SOURCE})
target_link_libraries(modify_wind gsl gslcblas m)

# inspect_wind
add_executable(inspect_wind ${PYTHON_SOURCE} ${INSPECT_WIND_SOURCE})
target_link_libraries(inspect_wind gsl gslcblas m)

# py_optd
add_executable(py_optd ${PYTHON_SOURCE} ${OPTICAL_DEPTH_SOURCE})
target_link_libraries(py_optd gsl gslcblas m)

# test
add_executable(py_unit_test ${TEST_SOURCE})
target_link_libraries(py_unit_test m gsl gslcblas cunit)
