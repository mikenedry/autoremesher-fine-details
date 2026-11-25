# Geometry regression tests

Build without Qt (requires installed TBB; Eigen is bundled):

```sh
cmake -S tests -B build/layout-tests -DCMAKE_BUILD_TYPE=Release
cmake --build build/layout-tests -j3
ctest --test-dir build/layout-tests --output-on-failure
```

Five executables cover the pipeline with synthetic geometry:

- `reference_surface_test`: immutable source provenance, disconnected parts,
  bounded refinement, diagonal balancing, affine grids, integer periods,
  source openings, thin tubes, separate vertex fans, and radial singularities.
- `surface_analysis_test`: curvature and scale covariance, directional sizing,
  feature transfer, protected curves and junctions, connected-sheet relaxation,
  and detection of missing surface area.
- `surface_relaxation_test`: local geometry relaxation and its numeric boundaries.
- `input_topology_test`: isolated negligible input fins, preserved open/closed
  geometry, outward branches and thin features, scale/order invariance, and original source IDs
  and boundary references after filtering. See the policy in
  [the headless documentation](../docs/HEADLESS.md#negligible-dangling-input-fins).
- `singularity_local_test`: exact agreement with the frozen published singularity
  implementation across 64 cases: open grids and closed tori, side lengths 12
  and 24, with/without feature edges, eight field seeds, and isolated vertices
  in the larger fixtures. Charges, before/after singularity counts, cancellation
  counts and every field vector must match. An instrumented allocator verifies
  fewer charge-vector-sized requests; timings are informational. The oracle's
  commit and hash are frozen in `reference/singularity_manifest.json` and checked
  at configuration time. Do not update it to match optimized production code.

The first executable also accepts an OBJ, a target quad count, and an output path.
It validates preparation provenance and basic output geometry before writing:

```sh
build/layout-tests/reference_surface_test input.obj 3000 output.obj
```

The root [headless build](../docs/HEADLESS.md) runs these same five programs against
one shared core target and adds the CLI contract test. That test exercises option
validation, source preservation, paths with spaces/non-ASCII characters, and an
actual procedural open-grid OBJ remesh with report/output consistency checks.
