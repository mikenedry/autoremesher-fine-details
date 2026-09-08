# Geometry regression tests

Build without Qt (requires installed TBB; Eigen is bundled):

```sh
cmake -S tests -B build/layout-tests -DCMAKE_BUILD_TYPE=Release
cmake --build build/layout-tests -j3
ctest --test-dir build/layout-tests --output-on-failure
```

Two executables cover the pipeline with synthetic geometry:

- `reference_surface_test`: immutable source provenance, disconnected parts,
  bounded refinement, diagonal balancing, affine grids, integer periods,
  source openings, thin tubes, separate vertex fans, and radial singularities.
- `surface_analysis_test`: curvature and scale covariance, directional sizing,
  feature transfer, protected curves and junctions, connected-sheet relaxation,
  and detection of missing surface area.

The first executable also accepts an OBJ, a target quad count, and an output path.
It validates preparation provenance and basic output geometry before writing:

```sh
build/layout-tests/reference_surface_test input.obj 3000 output.obj
```
