Unit Test Framework
###################

Unit tests are an excellent way to ensure that any code you write is robust and correct. To manage unit testing, SIROCCO
uses the `CUnit <https://gitlab.com/cunity/cunit>`_ test framework. Unit tests are run by using :code:`make check` in
either the root directory of SIROCCO, or in the source/test directory.

Installing CUnit
================

SIROCCO has been tested to work with CUnit (and CUnity) versions newer than 2.1-3. A recent version of CUnit is provided
in the :code:`$SIROCCO/software` directory and can be installed as a static library by using the Makefile in SIROCCO's
root directory. To build CUnit from source, you will need `CMake <https://cmake.org/>`_ installed, which is a modern
build system for C and C++ projects.

CUnit will be installed (as a static library) at the same time as GSL and SIROCCO during the first-time install, e.g.,

.. code:: bash

    $ [$SIROCCO] ./configure
    $ [$SIROCCO] make install

It is also possible to install only CUnit, using the same Makefile, if SIROCCO and GSL are already installed on your
system,

.. code:: bash

    $ [$SIROCCO] make cunit

If compilation of CUnit fails, it's more than likely that you could install a dynamic version of an older version of the
library from your system's package manager, e.g.

.. code:: bash

    # on macOS using homebrew
    $ brew install cunit

    # on Debian based Linux distributions
    $ sudo apt install libcunit1 libcunit1-doc libcunit1-dev

Running Tests
=============

To run the tests, navigate into one of three directories,

- :code:`$SIROCCO`
- :code:`$SIROCCO/source`
- :code:`$SIROCCO/tests`

Then run the command :code:`make check` which will compile and run the unit tests,

.. code:: bash

    $ [$SIROCCO/source] make check

    CUnit - A unit testing framework for C - Version 3.2.7-cunity
     http://cunit.sourceforge.net/

     Suite: Compton Processes
       Test: Klein-Nisina Formula ...passed
       Test: Compton Alpha - heating cross section  ...passed
       Test: Compton Beta - cooling cross section ...passed
       Test: Compton Formula ...passed
     Suite: Matrix Functions: GNU Science Library
       Test: Solve Matrix ...passed
       Test: Invert Matrix ...passed

     Run Summary       -      Run    Failed  Inactive   Skipped
          Suites       :        2         0         0         0
          Asserts      :       15         0       n/a       n/a
          Tests        :        6         0         0         0

     Elapsed Time: 0.009(s)

The important output is the "Run Summary", which lists the number of tests run and how many failed. To explain the
output a bit more, a suite is a collection of tests and a test is a function which evaluates the output of the function
being tested. If the output is deemed to be correct, the test passes. Otherwise, the function is counted as a failure.
The number of failed tests/suites is recorded in the "Failed" column of the table.

The asserts row in the table corresponds to the number of "checks" done, e.g. the number of function outputs checked for
correctness. In CUnit terminology, we assert that the output from a function should be something. If the output is that
something, then the test is a PASS otherwise it is marked as FAIL.

If a single test in a suite fails, the entire suite is marked as failed. In most cases, there are always more asserts
than suites and tests and there are always more, or an equal number of, tests than there are suites.

Writing Tests
=============

Creating a new test
-------------------

To create a test, we need to make a function which contains an assert statement from the CUnit library. An assert
statement is used to fail a test, so that if the condition in the assert statement is not true a failure is reported to
the CUnit test registry (more or that later). Test functions should not take any arguments and return an integer, which
is typically used to return an exit code which CUnit can use to determine is the test is successful or not it there are
no assert statements.

Assert statements come from the :code:`CUnit.h` header, with an exhaustive list of assertions available
`here <https://cunit.sourceforge.net/doc/writing_tests.html>`_. The code example below is a modified exert from the
one of matrix unit tests. In the function, test data is retrieved and compared to the output from :code:`solve_matrix`
using an assert which compares two floating point arrays to within a tolerance.

It should be noted that this assertion is not part of the standard CUnit assertions. It is possible to make a new
assertion by writing a macro (or function) which implements the base :code:`CU_assertImplementation` assert
implementation. If you need to create your own assertion, these should be kept in :code:`$SIROCCO/source/tests/assert.h`.

.. code:: c
    :caption: :code:`$SIROCCO/source/tests/tests/test_matrix.c`

    #include "assert.h"

    #include <CUnit/CUnit.h>

    int test_solve_matrix(void) {
      double *matrix_a;
      double *vector_b;
      double *vector_x;

      /* Get input data to `solve_matrix` and `vector_x` which is the "correct"
         answer we will use to compare to the output from `solve_matrix` */

      int vector_size;
      const int get_err =
        get_solve_matrix_test_data(..., &matrix_a, &vector_b, &vector_x, &vector_size);

      if (get_err) {  /* If we can't get the data, fail the test */
        CU_FAIL("Unable to load test data");  /* Assertion from CUnit.h */
      }

      /* Call `solve_matrix` with the input data from above */

      double *test_vector_x = malloc(vector_size * sizeof (double));
      const int matrix_err = solve_matrix(matrix_a, vector_b, vector_size, test_vector_x, -1);

      if (matrix_err) {  /* If there is some numerical error (or otherwise) fail the test */
        CU_FAIL("`solve_matrix` failed with error");
      }

      /* Use the following assertion to compare the value of the "correct" values (vector_x)
         against the output from `solve_matrix` (test_vector_x) */

      CU_ASSERT_DOUBLE_ARRAY_EQUAL_FATAL(test_vector_x, vector_x, vector_size, EPSILON);  /* Custom from assert.h */

      free(matrix_a);
      free(vector_b);
      free(vector_x);
      free(test_vector_x);

      return EXIT_SUCCESS;
    }


.. admonition:: Including :code:`python.h` in your tests

    If you need to access various structures or other things defined in :code:`python.h`, it is possible to include
    the header file in your test source code as in the example below (there are some data structures which depend
    on values defined in :code:`atomic.h`),

    .. code:: c

        #include "../../atomic.h"
        #include "../../python.h"

    In some situations this might complicate compilation of the unit test. In those cases, it could be better to
    re-define anything you need in the source file for the unit test.

Creating a test suite
---------------------

Unit tests belong in test suites and not by themselves. This means to create and run a unit test, we need a test suite
for that unit test to belong to. A test suite can be thought as a collection of tests, which are usually related. As an
example, there is a test suite for testing functions related to the Compton process and a test suite for matrix
functions.

The code exert below shows how to create a test suite and to add tests to the suite. The first step is to create a suite
to the CUnit test registry (the test registry is a global repository of test suites and associated tests) using
:code:`CU_add_suite`, which takes three arguments: the name of the suite, a function (pointer) to run when the suite
starts and a function to run after the suite has finished.

When a suite is added to the test registry, a pointer (:code:`CU_pSuite`) to the suite is returned from
:code:`CU_add_suite`. This pointer is used to add tests to the suite using :code:`CU_add_test` which takes three
arguments: a pointer to the suite to add the test to, the name of the test and the function (pointer) containing the
test. :code:`CU_add_test` returns a pointer to the test in the suite. If for whatever reason this fails, :code:`NULL` is
returned instead.

.. code:: c
    :caption: :code:`$SIROCCO/source/tests/tests/test_matrix.c`

    void create_matrix_test_suite(void) {
        /* Create a test suite - if suite can't be made, return error code */
        CU_pSuite suite = CU_add_suite(suite_name, matrix_suite_init, matrix_suite_teardown);
        if (suite == NULL) {
            CU_cleanup_registry();
            return CU_get_error();
        }

        /* Add some tests tests to suite - if one of them fails, return error code */
        if (CU_add_test(suite, "Solve Matrix", test_solve_matrix) == NULL) {
            CU_cleanup_registry();
            return CU_get_error();
        }
    }

The final two arguments for :code:`CU_add_suite` are used to initialise and clean up any additional data structures or
resources required to run the tests in the suite. In the matrix suite, for example, the cuSolver runtime is initialized
in `matrix_suite_init` and cleaned up in `matrix_suite_teardown`. An example of one of these functions, for the matrix
unit tests, is shown in the code exert below. These functions should not take any arguments and return an integer to
indicate if everything went OK or not.

.. code:: c
    :caption: :code:`$SIROCCO/source/tests/tests/test_matrix.c`



    int matrix_suite_init(void) {
        int error = EXIT_SUCCESS;

    #ifdef CUDA_ON  /* initialise cusolver */
        error = cusolver_create();
    #else  /* for GSL, we want to disable the default error handler */
        old_handler = gsl_set_error_handler_off();
    #endif

        return error;
    }

In the examples above, the code to create a suite and add tests is wrapped in a function
:code:`create_matrix_test_suite` with no arguments or return. All we need to do now to add those tests is to call that
function in the main function of the unit test framework, ensuring we do so after the test registry has been
initialized; this is done by the function  :code:`CU_initialize_registry`.

.. code:: c
    :caption: :code:`$SIROCCO/source/tests/unit_test_main.c`

    int main(int argc, char **argv) {
        /* Create the test registry */
        if (CU_initialize_registry() != CU_SUCCESS)   {
            return CU_get_error();
        }

        /* Add any test suites to the registry */
        create_matrix_test_suite();

        /* Set how verbose logging should be - CU_BRM_VERBOSE gets you the
           output shown in the running tests section */
        CU_basic_set_mode(CU_BRM_VERBOSE);

        /* Run the test suites */
        CU_basic_run_tests();

        /* Check how many tests failed */
        const int num_tests_failed = CU_get_number_of_tests_failed();

        /* Report on the number of tests failed, or if everything passed */
        if (num_tests_failed > 0) {
            printf("\033[1;31m%d test(s) failed\n\033[1;0m", num_tests_failed);  /* red text */
        } else {
            printf("\033[1;32mAll tests ran successfully\n\033[1;0m");  /* green text */
        }

        /* Clean up the CUnit registry */
        CU_cleanup_registry();

        return num_tests_failed;
    }

Directory and structure
-----------------------

Unit tests should be kept in logically named files within the unit test directory located at
:code:`$SIROCCO/source/tests/tests`. Any file in this directory should be added to the unit test Makefile, which is
located at :code:`$SIROCCO/source/tests/Makefile`, specifically to the :code:`TEST_SOURCES` variable which is a list of
all the source code required specifically for the unit test framework; this includes both the unit tests themselves and
any other code required to, e.g., build and control the test registry. Prototypes for wrapper functions for creating
test suites (which are called in the main function) should be placed in :code:`$SIROCCO/source/tests/tests/tests.h`
header file. Any data required for the tests should be kept in the data directory, :code:`$SIROCCO/source/tests/data`, in
appropriately organised directories as shown below.

.. code:: bash
    :caption: :code:`$SIROCCO/source/tests`

    $ tree $SIROCCO/source/tests

    ├── Makefile
    ├── assert.h
    ├── data
    │   └── matrix
    │       ├── inverse_macro
    │       │   ├── inverse.txt
    │       │   └── matrix.txt
    │       └── small_matrix
    │           ├── A.txt
    │           ├── b.txt
    │           └── x.txt
    ├── tests
    │   ├── test_matrix.c
    │   └── tests.h
    └── unit_test_main.c

We also need to include the SIROCCO source code we are testing in the :code:`SIROCCO_SOURCES` variable of the Makefile.
If there are any CUDA files required, these should be added to the :code:`CUDA_SOURCES` variable. In theory, we should
only need to include the files containing the code we are testing. But in practise, we choose to instead include all of
SIROCCO's source files (as it makes our lives easier) which increases compile time and the size of the final binary.
