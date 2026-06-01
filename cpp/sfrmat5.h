#pragma once

#include <array>
#include <memory>
#include <vector>

#include <Eigen/Dense>

namespace sfrmat5 {

enum class WindowFlag { Tukey = 0, Hamming = 1 };

template <typename T>
using Matrix = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;

template <typename T> struct SfrResult { // Outputs from slanted-edge SFR analysis.
    int status = 0;                      // 0 on success.
    Matrix<T> dat;                       // [frequency, mtf...] for each color/luminance.
    Matrix<T> e;                         // sampling efficiency (nval x ncol).
    T sfr50 = static_cast<T>(0);         // frequency where SFR = 50%.
    Matrix<T> fitme;                     // polynomial coefficients (+ misregistration if present).
    std::vector<T> esf;                  // last computed supersampled edge profile.
    int nbin = 4;                        // binning factor used.
    T del2 = static_cast<T>(0);          // sampling interval for ESF.
};

// ISO 12233 slanted-edge SFR analysis.
template <typename T> class SfrMat5 {
  public:
    SfrMat5();

    // RGB weights used to compute luminance (default: [0.213, 0.715, 0.072]).
    void set_weight(const std::array<T, 3>& weight);
    const std::array<T, 3>& weight() const;

    // Polynomial order for edge fit [1..5], default 5.
    void set_npol(int npol);
    int npol() const;

    // Window selection for edge/LSF processing.
    void set_wflag(WindowFlag wflag);
    WindowFlag wflag() const;

    // Sampling interval (mm or pixels/inch); if >1, treated as DPI.
    void set_del(T del);
    T del() const;

    // Compute SFR outputs for planar pixel data: channel 0 pixels, then channel 1, etc.
    SfrResult<T> compute(std::unique_ptr<std::vector<T>> pixels, int width, int height,
                         int channels) const;

  private:
    std::array<T, 3> weight_;              // RGB weights for luminance.
    int npol_ = 5;                         // polynomial edge-fit order.
    WindowFlag wflag_ = WindowFlag::Tukey; // window selection.
    T del_ = static_cast<T>(1);            // sampling interval (mm or DPI).
};

} // namespace sfrmat5
