# AGENTS

Use this guide to reproduce the C++ SFR computation results in this repo.

## Requirements

- `g++` with C++17 support
- `cmake`, `curl`, and `unzip` to build OpenCV locally
- Local OpenCV build in `third_party/opencv-install`

## Build and Run Tests

From repo root:

```bash
./scripts/build_opencv.sh
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
- The public matrix API uses `std::vector<std::vector<T>>`; OpenCV remains an internal dependency for least-squares polynomial fitting and DFT.
- `compute()` takes ownership of planar pixel data as `std::unique_ptr<std::vector<T>>` plus width, height, and channel count.
- Image loading/storage is outside the public SFR API; the test uses a local BMP helper to build planar pixel data.
