#pragma once

#include <array>
#include <vector>

namespace sfrmat5 {

struct Matrix {
    int rows = 0;
    int cols = 0;
    std::vector<double> data;

    Matrix() = default;
    Matrix(int r, int c, double value = 0.0);

    double &operator()(int r, int c);
    double operator()(int r, int c) const;
};

struct Image {
    int rows = 0;
    int cols = 0;
    int channels = 0;
    std::vector<double> data;

    Image() = default;
    Image(int r, int c, int ch, double value = 0.0);

    double &at(int r, int c, int ch);
    double at(int r, int c, int ch) const;
};

struct SfrResult {
    int status = 0;
    Matrix dat;    // [freq, mtf...]
    Matrix e;      // sampling efficiency (nval x ncol)
    double sfr50 = 0.0;
    Matrix fitme;  // polynomial coefficients (and misregistration if present)
    std::vector<double> esf; // last computed ESF
    int nbin = 4;
    double del2 = 0.0;
};

SfrResult compute_sfr(const Image &input,
                      double del = 1.0,
                      int npol = 5,
                      int wflag = 0,
                      const std::array<double, 3> &weight = {0.213, 0.715, 0.072});

}  // namespace sfrmat5
