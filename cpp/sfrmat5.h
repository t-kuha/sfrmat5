#pragma once

#include <array>
#include <vector>

namespace sfrmat5 {

template <typename T>
struct Matrix {
    int rows = 0;
    int cols = 0;
    std::vector<T> data;

    Matrix() = default;
    Matrix(int r, int c, T value = static_cast<T>(0));

    T &operator()(int r, int c);
    T operator()(int r, int c) const;
};

template <typename T>
struct Image {
    int rows = 0;
    int cols = 0;
    int channels = 0;
    std::vector<T> data;

    Image() = default;
    Image(int r, int c, int ch, T value = static_cast<T>(0));

    T &at(int r, int c, int ch);
    T at(int r, int c, int ch) const;
};

template <typename T>
struct SfrResult {
    int status = 0;
    Matrix<T> dat;    // [freq, mtf...]
    Matrix<T> e;      // sampling efficiency (nval x ncol)
    T sfr50 = static_cast<T>(0);
    Matrix<T> fitme;  // polynomial coefficients (and misregistration if present)
    std::vector<T> esf; // last computed ESF
    int nbin = 4;
    T del2 = static_cast<T>(0);
};

template <typename T>
class SfrMat5 {
public:
    SfrMat5();

    void set_weight(const std::array<T, 3> &weight);
    const std::array<T, 3> &weight() const;

    void set_npol(int npol);
    int npol() const;

    void set_wflag(int wflag);
    int wflag() const;

    void set_del(T del);
    T del() const;

    SfrResult<T> compute(const Image<T> &input) const;

private:
    std::array<T, 3> weight_;
    int npol_ = 5;
    int wflag_ = 0;
    T del_ = static_cast<T>(1);
};

}  // namespace sfrmat5
