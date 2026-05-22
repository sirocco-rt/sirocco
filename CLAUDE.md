# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

Sirocco (Simulating Ionization and Radiation in Outflows Created by Compact Objects) is a Monte Carlo radiative transfer code using the Sobolev approximation. It simulates winds and outflows in systems like cataclysmic variables, AGN, X-ray binaries, and young stellar objects. Formerly known as "Python" (renamed October 2024).

## Build Commands

**Environment setup** (required before building):
```bash
export SIROCCO=/path/to/Sirocco
```

**First-time full build** (compiles GSL, CUnit, then Sirocco):
```bash
./configure        # detects compilers (mpicc/gcc) and optional CUDA
make install       # builds GSL 2.6, CUnit 3.2.7, and all source
```

**Recompile source only** (after initial install):
```bash
cd source && make CC=mpicc sirocco       # main program
cd source && make CC=mpicc all           # all targets + indent check
cd source && make D sirocco              # debug/profiling build (-g -pg)
cd source && make CC=gcc sirocco         # compile without MPI
cd source && make INDENT=no sirocco      # skip auto-indentation
```

**Other build targets** (all from `source/`):
- `make swind` — spectral wind analysis
- `make windsave2table` — convert wind saves to ASCII tables
- `make windsave2fits` — convert wind saves to FITS (requires cfitsio)
- `make rad_hydro_files`, `make modify_wind`, `make inspect_wind` — wind utilities
- `make sirocco_optd` — optical depth calculations

**Unit tests** (requires CUnit + cmake):
```bash
make check                    # from top-level
cd source && make check       # from source directory
```

**Running Sirocco**:
```bash
sirocco parameter_file.pf                    # serial
mpirun -n 4 sirocco parameter_file.pf       # parallel with MPI
```

## Architecture

### Language and Compilation
- C (gnu99 standard), ~130 source files in `source/`
- Default compiler: `mpicc` (sets `-DMPI_ON`); falls back to gcc/clang
- Optional CUDA support via `--with-cuda` configure flag (matrix operations on GPU)
- Dependencies: GSL 2.6 (bundled in `software/`), CUnit 3.2.7 (bundled)
- Auto-indentation enforced on commit via `py_progs/run_indent.py` using GNU indent

### Key Source Files
- `source/sirocco.c` — main entry point, orchestrates ionization and spectrum cycles
- `source/sirocco.h` (~98KB) — central header with all major data structures (`PlasmaPtr`, `WindPtr`, `PhotPtr`, domain geometry structs)
- `source/atomic.h` (~31KB) — atomic data structures (ions, lines, levels, cross-sections)
- `source/templates.h` — auto-generated function prototypes (via `cproto`)
- `source/version.h` — auto-generated at build time from git hash

### Core Code Organization (in `source/`)
- **Wind models**: `define_wind.c`, `wind.c`, `wind2d.c`, `spherical.c`, `cylindrical.c`, `rtheta.c`, `sv.c` (Shlosman-Vitello), `knigge.c`, `homologous.c`, `corona.c`, `shell_wind.c`, `hydro_import.c`
- **Photon transport**: `trans_phot.c`, `photon2d.c`, `photon_gen.c`, `extract.c`, `paths.c`, `phot_util.c`
- **Radiation & spectra**: `radiation.c`, `spectra.c`, `bands.c`, `continuum.c`, `brem.c`, `compton.c`, `emission.c`
- **Ionization**: `ionization.c`, `direct_ion.c`, `saha.c`, `charge_exchange.c`, `partition.c`, `levels.c`, `recomb.c`
- **Macro-atom**: `matom.c`, `matom_diag.c`, `macro_gov.c`, `macro_gen_f.c`, `macro_accelerate.c`, `estimators_macro.c`
- **Line transfer**: `lines.c`, `resonance.c` (`resonate.c`), `dielectronic.c`
- **Disk**: `disk.c`, `disk_init.c`, `disk_photon_gen.c`
- **Setup/config**: `setup.c`, `setup_domains.c`, `setup_disk.c`, `setup_star_bh.c`, `setup_line_transfer.c`, `parse.c`, `rdpar.c`
- **MPI communication**: `communicate_plasma.c`, `communicate_wind.c`, `communicate_macro.c`, `communicate_spectra.c`, `para_update.c`
- **I/O**: `windsave.c`, `windsave2table_sub.c`, `xlog.c`, `diag.c`
- **Math utilities**: `recipes.c`, `random.c`, `cdf.c`, `vvector.c`, `matrix_cpu.c`, `matrix_gpu.cu`
- **Frame transformations**: `frame.c`

### Data Directories
- `xdata/` — atomic data files read at runtime (referenced via `Atomic_data` parameter)
- `xmod/` — model grids and spectra for disk/stellar models
- `examples/` — parameter files (`.pf`) organized by type: `basic/`, `extended/`, `regress/`, `gh-workflow/`

### Python Utilities (`py_progs/`)
Support scripts for data processing, visualization, and code maintenance:
- `run_indent.py` — enforces GNU indent code style on changed C files
- `MakeMacro.py` — generates macro-atom data from Chianti/Topbase databases
- `plot_spec.py`, `plot_wind.py`, `plot_tot.py` — visualization
- `hydro_2_python.py`, `import_1d.py`, `import_cyl.py`, `import_rtheta.py` — import external hydro models
- `balmer_decrement.py` — physics validation test
- `regression.py` — regression testing utilities

### Testing
- **Unit tests**: CUnit-based in `source/tests/` (test_matrix, test_compton, test_define_wind, test_run_mode, test_translate)
- **Integration tests**: GitHub Actions workflow (`.github/workflows/build.yml`) runs multiple example parameter files (CV, AGN, XRB, SN models) on pushes to dev/main
- **Regression tests**: Examples in `examples/regress/`

### Build System Notes
- `Makefile.in` is the top-level template (not auto-generated by autoconf — `configure` is a hand-written shell script)
- `source/Makefile` handles all C compilation, prototype generation, and indentation
- `source/tests/Makefile` includes the main source Makefile and links all Sirocco source for unit tests
- `get_models.c` cannot be included in `sirocco_source` list (prototype generation issue) but is added separately to object lists
- `make prototypes` regenerates `templates.h`, `log.h`, `atomic_proto.h`, `math_proto.h` via `cproto`

## Coding Conventions
- GNU indent style enforced automatically on changed files during build (unless `INDENT=no`)
- ANSI C with gnu99 standard
- MPI code guarded by `#ifdef MPI_ON` preprocessor directives
- CUDA code guarded by `#ifdef CUDA_ON`
- Logging via `xlog.c` functions (Log, Error, Debug)
