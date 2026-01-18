# AGENTS

Use this guide to reproduce the C++ SFR computation results in this repo.

## Requirements

- `g++` with C++17 support
- Vendored Eigen headers in `third_party/eigen-3.4.0`

## Build and Run Tests

From repo root:

```bash
./run_tests.sh
```

This compiles:

- `cpp/sfrmat5.cpp`
- `cpp/test_sfrmat5.cpp`

and executes the test against:

- `Example_Images/Test_edge1.bmp`

The test prints SFR50, sampling efficiency, and the first rows of SFR data.

## Notes

- The implementation is headless (no GUI/IO beyond the BMP loader in the test).
- To change the input image, edit `cpp/test_sfrmat5.cpp` and update the path.
- The C++ API is templated (`SfrMat5<T>`), with explicit instantiations for `float` and `double`.
- Use the instance accessors to set `weight`, `npol`, `wflag`, and `del` before calling `compute()`.
- `wflag` is a top-level enum (`WindowFlag`) with `Tukey` and `Hamming`.
- Eigen is used for least-squares polynomial fitting; FFT remains custom.
- Images are stored as per-channel Eigen matrices (`Image<T>::planes`); there is no `Image::at` accessor.
