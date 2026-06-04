#include "sfrmat5.h"

#define USE_MATH_DEFINES

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <vector>

#include <opencv2/core.hpp>

namespace sfrmat5 {

namespace {

template <typename T> struct Image {
    int rows = 0;
    int cols = 0;
    int channels = 0;
    std::vector<Matrix<T>> planes;

    /// Constructs an empty image container.
    Image() = default;

    /// Constructs an image with channel planes initialized to a constant value.
    Image(int r, int c, int ch, T value = static_cast<T>(0)) : rows(r), cols(c), channels(ch) {
        planes.reserve(ch);
        for (int i = 0; i < ch; ++i) {
            planes.emplace_back(r, c);
            planes.back().setTo(value);
        }
    }
};

template <typename T> Matrix<T> zeros(int rows, int cols) {
    return Matrix<T>(rows, cols, static_cast<T>(0));
}

double sum_region(const Matrix<double>& m, int row0, int row1, int col0, int col1) {
    double sum = 0.0;
    for (int r = row0; r < row1; ++r) {
        for (int c = col0; c < col1; ++c) {
            sum += m(r, c);
        }
    }
    return sum;
}

double mean_row(const Matrix<double>& m, int row) {
    return (m.cols == 0) ? 0.0 : sum_region(m, row, row + 1, 0, m.cols) / m.cols;
}

double mean_col(const Matrix<double>& m, int col) {
    return (m.rows == 0) ? 0.0 : sum_region(m, 0, m.rows, col, col + 1) / m.rows;
}

struct MeanStddev {
    double mean = 0.0;
    double stddev = 0.0;
};

/// Computes mean and sample standard deviation in one pass over the vector.
MeanStddev mean_stddev(const std::vector<double>& v) {
    if (v.empty()) {
        return {};
    }

    double sum = 0.0;
    double sum_squares = 0.0;
    for (double x : v) {
        sum += x;
        sum_squares += x * x;
    }

    MeanStddev stats;
    stats.mean = sum / static_cast<double>(v.size());
    if (v.size() > 1) {
        double sum_squared_deviations =
            sum_squares - static_cast<double>(v.size()) * stats.mean * stats.mean;
        stats.stddev = std::sqrt(std::max(0.0, sum_squared_deviations) /
                                 static_cast<double>(v.size() - 1));
    }
    return stats;
}

/// Returns the binomial coefficient n choose k as a double.
double nchoosek(int n, int k) {
    if (k < 0 || k > n) {
        return 0.0;
    }
    if (k == 0 || k == n) {
        return 1.0;
    }
    double res = 1.0;
    int kk = std::min(k, n - k);
    for (int i = 1; i <= kk; ++i) {
        res *= static_cast<double>(n - kk + i) / static_cast<double>(i);
    }
    return res;
}

/// Converts polynomial coefficients from scaled coordinates back to original x coordinates.
std::vector<double> polyfit_convert(const std::vector<double>& p2, const std::vector<double>& x) {
    int n = static_cast<int>(p2.size()) - 1;
    MeanStddev stats = mean_stddev(x);
    double m = stats.mean;
    double s = stats.stddev;
    if (s == 0.0) {
        s = 1.0;
    }
    std::vector<double> retval(p2.size(), 0.0);
    for (int i = 0; i <= n; ++i) {
        for (int j = 0; j <= i; ++j) {
            int idx = n - j;
            retval[idx] += p2[n - i] * nchoosek(i, j) * std::pow(-m, i - j) / std::pow(s, i);
        }
    }
    return retval;
}

/// Fits a polynomial after scaling x by its mean and sample standard deviation.
std::vector<double> polyfit_scaled(const std::vector<double>& x, const std::vector<double>& y,
                                   int degree) {
    if (x.size() != y.size()) {
        throw std::runtime_error("polyfit: x and y sizes differ");
    }
    int n = static_cast<int>(x.size());
    int m = degree + 1;
    MeanStddev stats = mean_stddev(x);
    double mx = stats.mean;
    double sx = stats.stddev;
    if (sx == 0.0) {
        sx = 1.0;
    }

    cv::Mat1d A(n, m);
    cv::Mat1d b(n, 1);
    for (int i = 0; i < n; ++i) {
        double z = (x[i] - mx) / sx;
        double value = 1.0;
        for (int p = m - 1; p >= 0; --p) {
            A(i, p) = value;
            value *= z;
        }
        b(i) = y[i];
    }

    cv::Mat1d p2;
    if (!cv::solve(A, b, p2, cv::DECOMP_SVD)) {
        throw std::runtime_error("polyfit: least-squares solve failed");
    }
    std::vector<double> coeffs(static_cast<size_t>(p2.rows), 0.0);
    for (int i = 0; i < p2.rows; ++i) {
        coeffs[i] = p2(i, 0);
    }
    return polyfit_convert(coeffs, x);
}

/// Evaluates a polynomial using Horner's method.
double polyval(const std::vector<double>& p, double x) {
    double y = 0.0;
    for (double coeff : p) {
        y = y * x + coeff;
    }
    return y;
}

/// Convolves two vectors and returns the centered output with the input length.
std::vector<double> conv_same(const std::vector<double>& x, const std::vector<double>& h) {
    int n = static_cast<int>(x.size());
    int m = static_cast<int>(h.size());
    std::vector<double> full(n + m - 1, 0.0);
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < m; ++j) {
            full[i + j] += x[i] * h[j];
        }
    }
    int start = m / 2;
    std::vector<double> same(n, 0.0);
    for (int i = 0; i < n; ++i) {
        same[i] = full[i + start];
    }
    return same;
}

/// Applies the finite-difference derivative filter to each matrix row.
Matrix<double> deriv1(const Matrix<double>& a, const std::vector<double>& fil) {
    Matrix<double> b = zeros<double>(a.rows, a.cols);
    for (int r = 0; r < a.rows; ++r) {
        std::vector<double> row(a.cols);
        for (int c = 0; c < a.cols; ++c) {
            row[c] = a(r, c);
        }
        std::vector<double> temp = conv_same(row, fil);
        for (int c = 0; c < a.cols; ++c) {
            b(r, c) = temp[c];
        }
        if (a.cols > 1) {
            b(r, 0) = b(r, 1);
            b(r, a.cols - 1) = b(r, a.cols - 2);
        }
    }
    return b;
}

/// Computes the intensity-weighted centroid location for a row profile.
double centroid(const std::vector<double>& x) {
    if (x.empty()) {
        return 0.0;
    }
    double sum = 0.0;
    double weighted_sum = 0.0;
    for (size_t i = 0; i < x.size(); ++i) {
        sum += x[i];
        weighted_sum += static_cast<double>(i + 1) * x[i];
    }
    if (sum == 0.0) {
        return 0.0;
    }
    return weighted_sum / sum;
}

/// Recenters a vector around the requested center location.
std::vector<double> cent(const std::vector<double>& a, double center) {
    int n = static_cast<int>(a.size());
    std::vector<double> b(n, 0.0);
    int mid = static_cast<int>(std::round((n + 1) / 2.0));
    int del = static_cast<int>(std::round(center - mid));
    if (del > 0) {
        for (int i = 0; i < n - del; ++i) {
            b[i] = a[i + del];
        }
    } else if (del < 0) {
        for (int i = -del; i < n; ++i) {
            b[i] = a[i + del];
        }
    } else {
        b = a;
    }
    return b;
}

/// Rotates an image 90 degrees counterclockwise.
Image<double> rotate90(const Image<double>& in) {
    Image<double> out(in.cols, in.rows, in.channels, 0.0);
    for (int ch = 0; ch < in.channels; ++ch) {
        cv::rotate(in.planes[ch], out.planes[ch], cv::ROTATE_90_COUNTERCLOCKWISE);
    }
    return out;
}

/// Rotates the image when the detected edge is closer to horizontal than vertical.
Image<double> rotatev2(const Image<double>& input) {
    Image<double> result = input;
    int nlin = input.rows;
    int npix = input.cols;
    int mm = (input.channels == 3 || input.channels == 4) ? 1 : 0;

    int nn = 3;
    int row_top = std::max(0, nn - 1);
    int row_bot = std::max(0, nlin - nn - 1);
    int col_left = std::max(0, nn - 1);
    int col_right = std::max(0, npix - nn - 1);

    const Matrix<double>& plane = input.planes[mm];
    double mean_bot = mean_row(plane, row_bot);
    double mean_top = mean_row(plane, row_top);
    double mean_right = mean_col(plane, col_right);
    double mean_left = mean_col(plane, col_left);

    double testv = std::abs(mean_bot - mean_top);
    double testh = std::abs(mean_right - mean_left);

    if (testv > testh) {
        result = rotate90(input);
    }
    return result;
}

/// Builds a Hamming window centered at the requested midpoint.
std::vector<double> ahamming(int n, double mid) {
    std::vector<double> data(n, 0.0);
    if (n == 0) {
        return data;
    }
    mid += 0.5;
    double wid1 = mid - 1.0;
    double wid2 = static_cast<double>(n) - mid;
    double wid = std::max(wid1, wid2);
    for (int i = 0; i < n; ++i) {
        double idx = static_cast<double>(i + 1);
        double arg = (idx - mid) * (M_PI / wid);
        data[i] = 0.54 + 0.46 * std::cos(arg);
    }
    return data;
}

/// Builds a symmetric Tukey window.
std::vector<double> tukey(int n, double alpha) {
    if (n == 1) {
        return {1.0};
    }
    if (alpha == 0.0) {
        return std::vector<double>(n, 1.0);
    }
    double m = (n - 1) / 2.0;
    int half = static_cast<int>(m);
    std::vector<double> wk(half + 1, 1.0);
    double thresh = alpha * m;
    for (int i = 0; i <= half; ++i) {
        double k = static_cast<double>(i);
        if (k <= thresh) {
            wk[i] = 0.5 * (1 + std::cos(M_PI * (k / (alpha * m) - 1)));
        }
    }
    std::vector<double> out(n, 0.0);
    for (int i = 0; i <= half; ++i) {
        out[i] = wk[i];
        out[n - 1 - i] = wk[i];
    }
    return out;
}

/// Builds a Tukey window shifted to the requested midpoint.
std::vector<double> tukey2(int n, double alpha, double mid) {
    if (n < 3) {
        return std::vector<double>(n, 1.0);
    }
    double m1 = n / 2.0;
    double m2 = mid;
    double m3 = n - mid;
    double mm = std::max(m2, m3);
    int n2 = static_cast<int>(std::round(2 * mm));
    std::vector<double> w = tukey(n2, alpha);
    if (mid >= m1) {
        w.resize(n, 0.0);
        return w;
    }
    int start = static_cast<int>(w.size()) - n;
    return std::vector<double>(w.begin() + start, w.begin() + start + n);
}

/// Computes correction factors for the derivative FIR frequency response.
std::vector<double> fir2fix(int n, int m) {
    std::vector<double> correct(n, 1.0);
    m = m - 1;
    for (int i = 1; i < n; ++i) {
        double num = M_PI * (i + 1) * m / (2.0 * (n + 1));
        double den = std::sin(M_PI * (i + 1) * m / (2.0 * (n + 1)));
        if (den == 0.0) {
            continue;
        }
        double val = std::abs(num / den);
        if (val > 10.0) {
            val = 10.0;
        }
        correct[i] = val;
    }
    return correct;
}

/// Fits the edge location data with a polynomial of the requested order.
std::vector<double> findedge2(const std::vector<double>& cent, int nlin, int nn) {
    std::vector<double> index(nlin, 0.0);
    for (int i = 0; i < nlin; ++i) {
        index[i] = static_cast<double>(i);
    }
    return polyfit_scaled(index, cent, nn);
}

struct ProjectResult {
    std::vector<double> point;
    int status = 0;
};

/// Projects a slanted edge image into a supersampled edge profile.
ProjectResult project2(const Matrix<double>& bb, const std::vector<double>& fitme, int fac) {
    int nlin = bb.rows;
    int npix = bb.cols;
    if (fac <= 0) {
        fac = 4;
    }

    double slope = fitme[fitme.size() - 2];
    slope = 1.0 / slope;
    int nn = static_cast<int>(std::floor(npix * fac));
    int offset = static_cast<int>(std::round(fac * (0 - (nlin - 1) / slope)));
    int del = std::abs(offset);
    if (offset > 0) {
        offset = 0;
    }
    int bwidth = nn + del + 150;
    std::vector<double> counts(bwidth, 0.0);
    std::vector<double> sums(bwidth, 0.0);

    std::vector<double> p2(nlin, 0.0);
    for (int m = 0; m < nlin; ++m) {
        double y = static_cast<double>(m);
        p2[m] = polyval(fitme, y) - fitme.back();
    }

    for (int n = 0; n < npix; ++n) {
        for (int m = 0; m < nlin; ++m) {
            double x = static_cast<double>(n);
            int ling = static_cast<int>(std::ceil((x - p2[m]) * fac)) + 1 - offset;
            if (ling < 1) {
                ling = 1;
            } else if (ling > bwidth) {
                ling = bwidth;
            }
            int idx = ling - 1;
            counts[idx] += 1.0;
            sums[idx] += bb(m, n);
        }
    }

    ProjectResult result;
    result.point.assign(nn, 0.0);
    int start = 1 + static_cast<int>(std::round(0.5 * del));

    for (int i = start; i < start + nn; ++i) {
        if (counts[i - 1] == 0.0) {
            result.status = 1;
            if (i == 1) {
                counts[i - 1] = counts[i];
                sums[i - 1] = sums[i];
            } else if (i == start + nn - 1) {
                counts[i - 1] = counts[i - 2];
                sums[i - 1] = sums[i - 2];
            } else {
                counts[i - 1] = (counts[i - 2] + counts[i]) / 2.0;
                sums[i - 1] = (sums[i - 2] + sums[i]) / 2.0;
            }
        }
    }

    for (int i = 0; i < nn; ++i) {
        int idx = i + start - 1;
        result.point[i] = sums[idx] / counts[idx];
    }
    return result;
}

/// Computes the discrete Fourier transform with OpenCV.
cv::Mat1d fft_magnitude(const std::vector<double>& x) {
    cv::Mat1d real(static_cast<int>(x.size()), 1);
    for (int i = 0; i < real.rows; ++i) {
        real(i, 0) = x[static_cast<size_t>(i)];
    }

    cv::Mat1d planes[] = {real, cv::Mat1d::zeros(real.rows, 1)};
    cv::Mat complex_signal;
    cv::merge(planes, 2, complex_signal);
    cv::dft(complex_signal, complex_signal);
    cv::split(complex_signal, planes);

    cv::Mat1d mag;
    cv::magnitude(planes[0], planes[1], mag);
    return mag;
}

/// Finds the spatial frequency where each SFR channel crosses the requested value.
std::vector<double> findfreq(const Matrix<double>& dat, double val, int imax, int fflag) {
    int nc = dat.cols - 1;
    std::vector<double> freqval(nc, 0.0);
    std::vector<double> sfrval(nc, 0.0);
    double maxf = dat(imax - 1, 0);
    std::vector<double> fil = {1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0};

    for (int c = 0; c < nc; ++c) {
        std::vector<double> col(dat.rows, 0.0);
        for (int r = 0; r < dat.rows; ++r) {
            col[r] = dat(r, c + 1);
        }
        if (fflag != 0) {
            std::vector<double> temp = conv_same(col, fil);
            for (int r = 1; r < dat.rows - 1; ++r) {
                col[r] = temp[r];
            }
        }

        int x = -1;
        for (int r = 0; r < imax; ++r) {
            if (col[r] - val < 0) {
                x = r - 1;
                break;
            }
        }
        double s = 0.0;
        double sval = 0.0;
        if (x <= 0) {
            s = maxf;
            sval = dat(imax - 1, c + 1);
        } else {
            sval = col[x];
            s = dat(x, 0);
            double y = col[x];
            double y2 = col[x + 1];
            double denom = (dat.rows > 1) ? dat(1, 0) : 0.0;
            double slope = (denom == 0.0) ? 0.0 : (y2 - y) / denom;
            double dely = col[x] - val;
            if (slope != 0.0) {
                s = s - dely / slope;
            }
            sval = sval - dely;
        }
        if (s > maxf) {
            s = maxf;
            sval = dat(imax - 1, c + 1);
        }
        freqval[c] = s;
        sfrval[c] = sval;
    }
    std::vector<double> out(2 * nc, 0.0);
    for (int c = 0; c < nc; ++c) {
        out[c] = freqval[c];
        out[nc + c] = sfrval[c];
    }
    return out;
}

/// Computes sampling efficiency percentages for requested SFR levels.
Matrix<double> sampeff(const Matrix<double>& dat, const std::vector<double>& val, double del,
                       int fflag) {
    if (dat.rows == 0 || dat.cols < 2) {
        return Matrix<double>();
    }
    double hs = 0.495 / del;
    int imax = dat.rows;
    int nindex = -1;
    for (int i = 0; i < dat.rows; ++i) {
        if (dat(i, 0) > hs) {
            nindex = i;
            break;
        }
    }
    if (nindex < 0) {
        Matrix<double> empty = zeros<double>(static_cast<int>(val.size()), dat.cols - 1);
        return empty;
    }

    int nc = dat.cols - 1;
    Matrix<double> eff = zeros<double>(static_cast<int>(val.size()), nc);
    for (size_t v = 0; v < val.size(); ++v) {
        std::vector<double> freq_sfr = findfreq(dat, val[v], imax, fflag);
        for (int c = 0; c < nc; ++c) {
            double freq = std::clamp(freq_sfr[c], 0.0, hs);
            eff(static_cast<int>(v), c) = std::min(std::round(100.0 * freq / hs), 100.0);
        }
    }
    return eff;
}

/// Computes coefficient of determination and root mean square error.
void rsquare(const std::vector<double>& y, const std::vector<double>& f, double& r2, double& rmse) {
    if (y.size() != f.size() || y.empty()) {
        r2 = 0.0;
        rmse = 0.0;
        return;
    }
    std::vector<double> yy;
    std::vector<double> ff;
    for (size_t i = 0; i < y.size(); ++i) {
        if (!std::isnan(y[i]) && !std::isnan(f[i])) {
            yy.push_back(y[i]);
            ff.push_back(f[i]);
        }
    }
    double mean_y = mean_stddev(yy).mean;
    double ss_res = 0.0;
    double ss_tot = 0.0;
    for (size_t i = 0; i < yy.size(); ++i) {
        double diff = yy[i] - ff[i];
        ss_res += diff * diff;
        double dt = yy[i] - mean_y;
        ss_tot += dt * dt;
    }
    r2 = (ss_tot == 0.0) ? 0.0 : std::max(0.0, 1.0 - ss_res / ss_tot);
    rmse = std::sqrt(ss_res / static_cast<double>(yy.size()));
}

/// Runs the double-precision SFR pipeline and returns raw double outputs.
SfrResult<double> compute_sfr_double(const Image<double>& input, double del, int npol,
                                     WindowFlag wflag, const std::array<double, 3>& weight) {
    if (input.rows == 0 || input.cols == 0) {
        throw std::runtime_error("Empty input image");
    }
    const int nbin = 4;
    const double alpha = 1.0;
    npol = std::min(npol, 5);

    Image<double> a = input;
    if (a.channels == 3) {
        Image<double> out(a.rows, a.cols, 4, 0.0);
        out.planes[0] = a.planes[0];
        out.planes[1] = a.planes[1];
        out.planes[2] = a.planes[2];
        out.planes[3] = weight[0] * a.planes[0] + weight[1] * a.planes[1] + weight[2] * a.planes[2];
        a = out;
    }

    if (del > 1.0) {
        del = 25.4 / del;
    }

    a = rotatev2(a);
    int nlin = a.rows;
    int npix = a.cols;
    int ncol = a.channels;

    int left_cols = std::min(5, npix);
    int right_cols = std::min(6, npix);
    const Matrix<double>& plane0 = a.planes[0];
    double tleft = sum_region(plane0, 0, plane0.rows, 0, left_cols);
    double tright = sum_region(plane0, 0, plane0.rows, plane0.cols - right_cols, plane0.cols);

    std::vector<double> fil1 = {0.5, -0.5};
    std::vector<double> fil2 = {0.5, 0.0, -0.5};
    if (tleft > tright) {
        fil1 = {-0.5, 0.5};
        fil2 = {-0.5, 0.0, 0.5};
    }

    std::vector<double> win1;
    if (wflag == WindowFlag::Hamming) {
        win1 = ahamming(npix, (npix + 1) / 2.0);
    } else {
        win1 = tukey2(npix, alpha, (npix + 1) / 2.0);
        for (double& v : win1) {
            v = v * 0.95 + 0.05;
        }
    }

    std::vector<std::vector<double>> loc(ncol, std::vector<double>(nlin, 0.0));
    std::vector<std::vector<double>> fitme(ncol);
    std::vector<std::vector<double>> fitme1(ncol);

    for (int color = 0; color < ncol; ++color) {
        Matrix<double> plane = a.planes[color];
        Matrix<double> deriv = deriv1(plane, fil1);
        for (int n = 0; n < nlin; ++n) {
            std::vector<double> row(npix, 0.0);
            for (int i = 0; i < npix; ++i) {
                row[i] = deriv(n, i) * win1[i];
            }
            loc[color][n] = centroid(row) - 0.5;
        }
        fitme[color] = findedge2(loc[color], nlin, npol);

        for (int n = 0; n < nlin; ++n) {
            double place = polyval(fitme[color], static_cast<double>(n));
            std::vector<double> win2 =
                (wflag == WindowFlag::Hamming) ? ahamming(npix, place) : tukey2(npix, alpha, place);
            if (wflag == WindowFlag::Tukey) {
                for (double& v : win2) {
                    v = v * 0.95 + 0.05;
                }
            }
            std::vector<double> row(npix, 0.0);
            for (int i = 0; i < npix; ++i) {
                row[i] = deriv(n, i) * win2[i];
            }
            loc[color][n] = centroid(row) - 0.5;
        }

        fitme[color] = findedge2(loc[color], nlin, npol);
        fitme1[color] = findedge2(loc[color], nlin, 1);

        if (npol > 3) {
            std::vector<double> x(nlin, 0.0);
            std::vector<double> y(nlin, 0.0);
            for (int i = 0; i < nlin; ++i) {
                x[i] = i;
                y[i] = polyval(fitme[color], x[i]);
            }
            double r2 = 0.0;
            double rmse = 0.0;
            rsquare(y, loc[color], r2, rmse);
        }
    }

    std::vector<double> midloc(ncol, 0.0);
    std::vector<double> misreg(ncol, 0.0);
    for (int i = 0; i < ncol; ++i) {
        midloc[i] = polyval(fitme[i], (nlin - 1) / 2.0);
    }
    if (ncol > 2) {
        for (int i = 0; i < ncol; ++i) {
            misreg[i] = midloc[i] - midloc[1];
        }
    }

    double slope_ref = fitme1.back()[fitme1.back().size() - 2];
    int nlin1 = nlin;
    if (std::abs(slope_ref) > std::numeric_limits<double>::epsilon()) {
        nlin1 = static_cast<int>(
            std::round(std::floor(nlin * std::abs(slope_ref)) / std::abs(slope_ref)));
    }
    nlin1 = std::max(1, std::min(nlin1, nlin));
    if (nlin1 < nlin) {
        Image<double> cropped(nlin1, npix, ncol, 0.0);
        for (int ch = 0; ch < ncol; ++ch) {
            cropped.planes[ch] = a.planes[ch].rowRange(0, nlin1).clone();
        }
        a = cropped;
        nlin = a.rows;
    }

    double vslope = -slope_ref;
    double delfac = std::cos(std::atan(vslope));
    double delimage = del;
    del = del * delfac;
    double del2 = del / nbin;
    if (ncol > 2) {
        for (double& m : misreg) {
            m *= delfac;
        }
    }

    int nn = static_cast<int>(std::ceil(npix * nbin));
    int nn2 = static_cast<int>(std::floor(nn / 2.0)) + 1;
    std::vector<double> dcorr = fir2fix(nn2, 3);
    int freqlim = (nbin == 1) ? 2 : 1;
    int nn2out = static_cast<int>(std::round(nn2 * freqlim / 2.0));

    Matrix<double> mtf = zeros<double>(nn, ncol);
    std::vector<double> esf_last;

    for (int color = 0; color < ncol; ++color) {
        Matrix<double> plane = a.planes[color];
        ProjectResult proj = project2(plane, fitme[color], nbin);
        esf_last = proj.point;
        Matrix<double> esf_mat = zeros<double>(1, static_cast<int>(esf_last.size()));
        for (int i = 0; i < static_cast<int>(esf_last.size()); ++i) {
            esf_mat(0, i) = esf_last[i];
        }
        Matrix<double> deriv = deriv1(esf_mat, fil2);
        std::vector<double> c(nn, 0.0);
        for (int i = 0; i < nn; ++i) {
            c[i] = deriv(0, i);
        }
        if (!c.empty()) {
            if (c.front() == 0.0 && c.size() > 1) {
                c.front() = c[1];
            } else if (c.back() == 0.0 && c.size() > 1) {
                c.back() = c[c.size() - 2];
            }
        }

        double max_val = *std::max_element(c.begin(), c.end());
        double mm = 0.0;
        int count = 0;
        for (int i = 0; i < nn; ++i) {
            if (c[i] == max_val) {
                mm += i + 1;
                count++;
            }
        }
        if (count > 0) {
            mm /= count;
        }
        c = cent(c, mm);
        double center = nn / 2.0;
        std::vector<double> win =
            (wflag == WindowFlag::Hamming) ? ahamming(nn, center) : tukey2(nn, alpha, center);
        for (int i = 0; i < nn; ++i) {
            c[i] *= win[i];
        }

        cv::Mat1d fx = fft_magnitude(c);
        double dc0 = fx(0, 0);
        for (int i = 0; i < nn2; ++i) {
            double val = (dc0 == 0.0) ? 0.0 : fx(i, 0) / dc0;
            val *= dcorr[i];
            mtf(i, color) = val;
        }
    }

    std::vector<double> freq(nn, 0.0);
    for (int i = 0; i < nn; ++i) {
        freq[i] = static_cast<double>(i) / (del2 * nn);
    }
    Matrix<double> dat = zeros<double>(nn2out, ncol + 1);
    for (int r = 0; r < nn2out; ++r) {
        dat(r, 0) = freq[r];
        for (int cidx = 0; cidx < ncol; ++cidx) {
            dat(r, cidx + 1) = mtf(r, cidx);
        }
    }

    int fit_cols = static_cast<int>(fitme[0].size());
    Matrix<double> fitout = zeros<double>(ncol, (ncol > 2) ? fit_cols + 1 : fit_cols);
    for (int r = 0; r < ncol; ++r) {
        for (int c = 0; c < fit_cols; ++c) {
            fitout(r, c) = fitme[r][c];
        }
        if (ncol > 2) {
            fitout(r, fit_cols) = misreg[r];
        }
    }

    std::vector<double> val = {0.1, 0.5};
    Matrix<double> eff = sampeff(dat, val, delimage, 0);
    std::vector<double> freq_sfr = findfreq(dat, 0.5, dat.rows, 0);
    double sfr50 = freq_sfr.empty() ? 0.0 : freq_sfr[0];

    SfrResult<double> result;
    result.status = 0;
    result.dat = dat;
    result.e = eff;
    result.sfr50 = sfr50;
    result.fitme = fitout;
    result.esf = esf_last;
    result.nbin = nbin;
    result.del2 = del2;
    return result;
}

/// Converts planar input pixels to the internal double-precision image representation.
template <typename T>
Image<double> to_double_image(const std::vector<T>& pixels, int width, int height, int channels) {
    if (width <= 0 || height <= 0 || channels <= 0) {
        throw std::invalid_argument("SfrMat5::compute requires positive dimensions");
    }

    const auto expected_size =
        static_cast<size_t>(width) * static_cast<size_t>(height) * static_cast<size_t>(channels);
    if (pixels.size() != expected_size) {
        throw std::invalid_argument("SfrMat5::compute pixel data size does not match dimensions");
    }

    Image<double> out(height, width, channels, 0.0);
    const size_t plane_size = static_cast<size_t>(width) * static_cast<size_t>(height);
    for (int ch = 0; ch < channels; ++ch) {
        const size_t channel_offset = static_cast<size_t>(ch) * plane_size;
        for (int row = 0; row < height; ++row) {
            const size_t row_offset = channel_offset + static_cast<size_t>(row) * width;
            for (int col = 0; col < width; ++col) {
                out.planes[ch](row, col) = static_cast<double>(pixels[row_offset + col]);
            }
        }
    }
    return out;
}

/// Casts a double-precision matrix to the requested scalar type.
template <typename T> Matrix<T> cast_matrix(const Matrix<double>& input) {
    Matrix<T> out(input.rows, input.cols);
    for (int r = 0; r < input.rows; ++r) {
        for (int c = 0; c < input.cols; ++c) {
            out(r, c) = static_cast<T>(input(r, c));
        }
    }
    return out;
}

/// Casts a double-precision vector to the requested scalar type.
template <typename T> std::vector<T> cast_vector(const std::vector<double>& input) {
    std::vector<T> out(input.size(), static_cast<T>(0));
    for (size_t i = 0; i < input.size(); ++i) {
        out[i] = static_cast<T>(input[i]);
    }
    return out;
}

/// Casts a double-precision SFR result to the requested scalar type.
template <typename T> SfrResult<T> cast_result(const SfrResult<double>& input) {
    SfrResult<T> out;
    out.status = input.status;
    out.dat = cast_matrix<T>(input.dat);
    out.e = cast_matrix<T>(input.e);
    out.sfr50 = static_cast<T>(input.sfr50);
    out.fitme = cast_matrix<T>(input.fitme);
    out.esf = cast_vector<T>(input.esf);
    out.nbin = input.nbin;
    out.del2 = static_cast<T>(input.del2);
    return out;
}

} // namespace

/// Constructs an analyzer with default luminance weights, polynomial order, window, and sampling.
template <typename T>
SfrMat5<T>::SfrMat5()
    : weight_{static_cast<T>(0.213), static_cast<T>(0.715), static_cast<T>(0.072)}, npol_(5),
      wflag_(WindowFlag::Tukey), del_(static_cast<T>(1)) {}

/// Sets RGB weights used to compute luminance.
template <typename T> void SfrMat5<T>::set_weight(const std::array<T, 3>& weight) {
    weight_ = weight;
}

/// Returns RGB weights used to compute luminance.
template <typename T> const std::array<T, 3>& SfrMat5<T>::weight() const {
    return weight_;
}

/// Sets polynomial order for edge fit, clamped to [1, 5].
template <typename T> void SfrMat5<T>::set_npol(int npol) {
    // Ensure npol is within valid range [1, 5]
    npol_ = std::clamp(npol, 1, 5);
}

/// Returns polynomial order for edge fit.
template <typename T> int SfrMat5<T>::npol() const {
    return npol_;
}

/// Sets window selection for edge and LSF processing.
template <typename T> void SfrMat5<T>::set_wflag(WindowFlag wflag) {
    wflag_ = wflag;
}

/// Returns window selection for edge and LSF processing.
template <typename T> WindowFlag SfrMat5<T>::wflag() const {
    return wflag_;
}

/// Sets sampling interval in millimeters, or DPI when greater than 1.
template <typename T> void SfrMat5<T>::set_del(T del) {
    del_ = del;
}

/// Returns sampling interval setting.
template <typename T> T SfrMat5<T>::del() const {
    return del_;
}

/// Computes SFR outputs from planar pixels.
template <typename T>
SfrResult<T> SfrMat5<T>::compute(std::unique_ptr<std::vector<T>> pixels, int width, int height,
                                 int channels) const {
    if (!pixels) {
        throw std::invalid_argument("SfrMat5::compute requires non-null pixel data");
    }
    Image<double> img = to_double_image(*pixels, width, height, channels);
    std::array<double, 3> w = {static_cast<double>(weight_[0]), static_cast<double>(weight_[1]),
                               static_cast<double>(weight_[2])};
    SfrResult<double> res = compute_sfr_double(img, static_cast<double>(del_), npol_, wflag_, w);
    return cast_result<T>(res);
}

template class SfrMat5<float>;
template class SfrMat5<double>;

template struct SfrResult<float>;
template struct SfrResult<double>;

} // namespace sfrmat5
