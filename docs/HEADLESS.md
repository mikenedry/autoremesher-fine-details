# Standalone library and OBJ driver

The root CMake build provides `autoremesher_core` (`AutoRemesher::core` for CMake
consumers) and the Qt-free OBJ driver `AutoRemesherEngine`.

## Build and test

Requirements: a C++14 compiler, CMake 3.20 or newer, and an installed TBB CMake
package. Eigen, isotropicremesher, meshoptimizer and tinyobjloader come from this
source tree. macOS uses Accelerate and TBB's allocator targets, retaining the
existing headless build's compiler and solver settings. Tests also need Python 3
(standard library only for the CLI test).

```sh
cmake -S . -B build/headless -DCMAKE_BUILD_TYPE=Release
cmake --build build/headless --parallel 2
ctest --test-dir build/headless --output-on-failure
build/headless/AutoRemesherEngine --version
```

For a multi-configuration generator, pass `--config Release` to the build,
`-C Release` to CTest, and use the executable under `build/headless/Release`.
`AUTOREMESHER_BUILD_CLI=OFF` builds just the library and selected tests;
`AUTOREMESHER_BUILD_TESTS=OFF` removes the test/Python dependency. A parent CMake
project can use `add_subdirectory` and link `AutoRemesher::core`. Installation and
relocation of this SDK are not yet qualified.

The version identifies the configured Git revision, adding `+dirty` for local
changes. For an extracted source archive, set `AUTOREMESHER_SOURCE_REVISION` to its
recorded commit; otherwise the revision is explicitly `unknown`. Reconfigure
after changing commits. Release packaging must verify the source and artifact
hashes; this version string alone is not a provenance guarantee.

## OBJ invocation and compatibility

```sh
build/headless/AutoRemesherEngine --input input.obj --output output.obj \
  --quads 5000 --no-cloth-folds --report result.json --quiet
```

The driver implements version 2 of the OBJ command-line protocol. Its
cloth-fold guidance defaults to **enabled**, unlike the engine API and Qt
application defaults; `--no-cloth-folds` is the explicit opt-out. This driver uses
`--quads`, while the existing Qt application's batch mode uses `--target-quads`.

The driver triangulates OBJ input and writes indexed polygon output. The JSON
report records actual vertex, face and true-quad counts, cloth-fold selection,
source revision and elapsed seconds. The target is a request, not an exact-count
guarantee. The legacy progress lines and final JSON share stdout; `--quiet`
suppresses the driver's progress callback. Diagnostics and the engine phase
report use stderr. Exit 0 currently means a nonempty result passed the driver's
basic output checks; exit 2 reports a caught option, input, engine or IO error.
The driver does not check completeness per input component or implement cooperative cancellation.

The marker `local-singularity-1` identifies the imported local-charge
optimization. It replaces repeated whole-mesh charge scans during pair
cancellation with queries for the affected vertices. Pair search, relaxation,
acceptance and rollback rules are unchanged. The regression compiles a frozen,
hash-checked copy of the published implementation as a separate oracle and
compares charges, cancellation counts and every field vector exactly. It also
checks reduced charge-vector-sized allocation requests; elapsed times are
informational.

## Source and licensing

The OBJ driver, the local-charge changes in
`src/AutoRemesher/singularitysimplifier.cpp` and the corresponding
`tests/singularity_local_test.cpp` carry the
[GPL-3.0-or-later license](../LICENSES/GPL-3.0-or-later.txt).
The original source retains its [MIT notice](../LICENSE). The frozen test
reference is the original MIT source at commit
`fef338d27aeb9c5f0eca12b305538ed670c59776`. Its hash is recorded in
`tests/reference/singularity_manifest.json`. Third-party components retain
their own notices.

## Negligible dangling input fins

Before island separation, the core recognizes a narrow class of non-manifold
input artifacts: a single triangle attached to an otherwise continuous sheet.
It qualifies only when two edges are free boundaries, the attachment edge has
exactly three incident faces, and its free vertex belongs to no other face.
The other two faces must have consistent opposite edge traversal and supporting
normals within five degrees, with their third vertices on opposite sides.

The fin tip must project strictly inside the attachment segment. Its distance
from that segment must be at most 0.5% of the segment length and at most 0.2%
of either supporting triangle's altitude. These local, scale-independent bounds
limit the geometry approximation; small area alone never permits removal.
For nonzero-height fins, the tip must also lie inward between the two support
planes: on each plane's same side as the opposite supporting vertex. Outward
branches are preserved. Flat or numerically unresolved supporting bends have no
reliable inward side and are preserved too. This test is independent of global
winding; exactly zero-height degenerate fins retain the preceding bounds.
This policy cannot infer arbitrary artist intent: a deliberately authored internal fin
matching all these negligible-branch criteria is treated as redundant. Larger
branches, connected tips, creases, ordinary open boundaries and manifold thin
triangles are preserved. Invalid or ambiguous input is not repaired by this rule.

Only qualifying faces are omitted from a private working input. Retained vertex
positions, face order and original source IDs are preserved. Sizing, island
separation, surface analysis and boundary references all use the same corrected
input; the original constructor input remains unchanged. The core exposes
`removedInputFinTriangles()` as original face IDs and prints a removal count.
No output hole filling or change to the field solver/extractor is involved.

`input_topology_test` covers the three-face regression, original-ID/reference
consistency, ordinary closed/open surfaces, authored openings, thin triangles,
outward/substantial/connected branches, flat or inconsistent supports, geometric bounds and
invariance under face ordering, winding, rotation and scales from 1e-9 to 1e9.

## Unicode CLI paths

On Windows, the driver receives UTF-16 arguments through `wmain` and opens
input OBJ, output OBJ and report paths with the Windows C++ runtime's wide
filename stream overloads. The CLI remains C++14; MinGW builds link the CLI
with `-municode`. Non-Windows paths keep their native byte representation.
Option text and diagnostics use UTF-8 on Windows.

The input stream goes through the same tinyobjloader parser, triangulation and
vertex-color defaults. Its `MaterialFileReader` remains enabled with the prior
empty base directory (material paths resolve from the working directory).
Encoding and lookup of filenames embedded in OBJ material references are
unchanged; this fix concerns the three command-line filesystem paths.

`cli_paths` exercises ASCII, accented, CJK and supplementary Unicode names,
spaces, empty-input validation, and synthetic OBJ input/output/report round trips.
It can also be run independently with `python tests/cli_paths_test.py /path/to/AutoRemesherEngine`.
