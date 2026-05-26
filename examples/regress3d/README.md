# 3D Regression Tests (x3d branch)

This directory contains regression tests for the CYLIND3D (true 3D cylindrical)
geometry introduced on the x3d branch.

Models here are intended to verify that the CYLIND3D wind model initialises,
transports photons, and writes wind-saves correctly.  Like the standard
regression tests they need not be converged — the goal is to catch gross
regressions in 3D-specific code paths.
