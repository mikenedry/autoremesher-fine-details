# Geometry regression tests

Build without Qt using bundled Eigen and installed TBB:

```sh
cmake -S tests -B build/layout-tests -DCMAKE_BUILD_TYPE=Release
cmake --build build/layout-tests -j3
ctest --test-dir build/layout-tests --output-on-failure
```

The preparation test also validates and remeshes an external OBJ:

```sh
build/layout-tests/reference_surface_test input.obj 3000 output.obj
```
