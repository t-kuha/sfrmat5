#pragma once

#include <array>
#include <memory>
#include <vector>

namespace sfrmat5 {

enum class WindowFlag { Tukey = 0, Hamming = 1 };

/// Row-major matrix used by the public SFR API.
template <typename T> using Matrix = std::vector<std::vector<T>>;

/// Stores outputs from slanted-edge SFR analysis.
template <typename T> struct SfrResult {
    int status = 0;              // 0 on success.
    Matrix<T> dat;               // [frequency, mtf...] for each color/luminance.
    Matrix<T> e;                 // sampling efficiency (nval x ncol).
    T sfr50 = static_cast<T>(0); // frequency where SFR = 50%.
    Matrix<T> fitme;             // polynomial coefficients (+ misregistration if present).
    std::vector<T> esf;          // last computed supersampled edge profile.
    int nbin = 4;                // binning factor used.
    T del2 = static_cast<T>(0);  // sampling interval for ESF.
};

/// Performs ISO 12233 slanted-edge SFR analysis on planar pixel data.
template <typename T> class SfrMat5 {
  public:
    /// Constructs an analyzer with MATLAB-compatible defaults.
    SfrMat5();

    /// Sets RGB weights used to compute luminance.
    void set_weight(const std::array<T, 3>& weight);

    /// Returns RGB weights used to compute luminance.
    const std::array<T, 3>& weight() const;

    /// Sets polynomial order for edge fit, clamped to [1, 5].
    void set_npol(int npol);

    /// Returns polynomial order for edge fit.
    int npol() const;

    /// Sets window selection for edge and LSF processing.
    void set_wflag(WindowFlag wflag);

    /// Returns window selection for edge and LSF processing.
    WindowFlag wflag() const;

    /// Sets sampling interval in millimeters, or DPI when greater than 1.
    void set_del(T del);

    /// Returns sampling interval setting.
    T del() const;

    /// Computes SFR outputs from planar pixels: all channel 0 pixels, then channel 1, etc.
    SfrResult<T> compute(std::unique_ptr<std::vector<T>> pixels, int width, int height,
                         int channels) const;

  private:
    std::array<T, 3> weight_;              // RGB weights for luminance.
    int npol_ = 5;                         // polynomial edge-fit order.
    WindowFlag wflag_ = WindowFlag::Tukey; // window selection.
    T del_ = static_cast<T>(1);            // sampling interval (mm or DPI).
};

} // namespace sfrmat5
