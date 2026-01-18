#include "sfrmat5.h"

#define USE_MATH_DEFINES

#include <algorithm>
#include <cmath>
#include <complex>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <vector>

namespace sfrmat5 {

Matrix::Matrix(int r, int c, double value) : rows(r), cols(c), data(r * c, value) {}
double &Matrix::operator()(int r, int c) { return data[r * cols + c]; }
double Matrix::operator()(int r, int c) const { return data[r * cols + c]; }

Image::Image(int r, int c, int ch, double value)
    : rows(r), cols(c), channels(ch), data(r * c * ch, value) {}
double &Image::at(int r, int c, int ch) { return data[(ch * rows + r) * cols + c]; }
double Image::at(int r, int c, int ch) const { return data[(ch * rows + r) * cols + c]; }

namespace {

double clip_value(double in, double low, double high) {
    return std::min(std::max(in, low), high);
}

std::vector<double> clip_vector(const std::vector<double> &in, double low, double high) {
    std::vector<double> out = in;
    for (double &v : out) {
        v = clip_value(v, low, high);
    }
    return out;
}

double mean(const std::vector<double> &v) {
    if (v.empty()) {
        return 0.0;
    }
    double sum = std::accumulate(v.begin(), v.end(), 0.0);
    return sum / static_cast<double>(v.size());
}

double stddev(const std::vector<double> &v, double m) {
    if (v.size() < 2) {
        return 0.0;
    }
    double acc = 0.0;
    for (double x : v) {
        double d = x - m;
        acc += d * d;
    }
    return std::sqrt(acc / static_cast<double>(v.size() - 1));
}

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

std::vector<double> polyfit_convert(const std::vector<double> &p2, const std::vector<double> &x) {
    int n = static_cast<int>(p2.size()) - 1;
    double m = mean(x);
    double s = stddev(x, m);
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

std::vector<double> solve_linear(std::vector<std::vector<double>> a, std::vector<double> b) {
    int n = static_cast<int>(b.size());
    for (int i = 0; i < n; ++i) {
        int pivot = i;
        double max_val = std::abs(a[i][i]);
        for (int r = i + 1; r < n; ++r) {
            if (std::abs(a[r][i]) > max_val) {
                max_val = std::abs(a[r][i]);
                pivot = r;
            }
        }
        if (max_val == 0.0) {
            throw std::runtime_error("Singular matrix in polyfit");
        }
        if (pivot != i) {
            std::swap(a[pivot], a[i]);
            std::swap(b[pivot], b[i]);
        }
        double diag = a[i][i];
        for (int c = i; c < n; ++c) {
            a[i][c] /= diag;
        }
        b[i] /= diag;
        for (int r = 0; r < n; ++r) {
            if (r == i) {
                continue;
            }
            double factor = a[r][i];
            for (int c = i; c < n; ++c) {
                a[r][c] -= factor * a[i][c];
            }
            b[r] -= factor * b[i];
        }
    }
    return b;
}

std::vector<double> polyfit_scaled(const std::vector<double> &x,
                                   const std::vector<double> &y,
                                   int degree) {
    if (x.size() != y.size()) {
        throw std::runtime_error("polyfit: x and y sizes differ");
    }
    int n = static_cast<int>(x.size());
    int m = degree + 1;
    double mx = mean(x);
    double sx = stddev(x, mx);
    if (sx == 0.0) {
        sx = 1.0;
    }

    std::vector<std::vector<double>> ata(m, std::vector<double>(m, 0.0));
    std::vector<double> aty(m, 0.0);

    for (int i = 0; i < n; ++i) {
        double z = (x[i] - mx) / sx;
        std::vector<double> v(m, 1.0);
        for (int p = 1; p < m; ++p) {
            v[p] = v[p - 1] * z;
        }
        for (int row = 0; row < m; ++row) {
            for (int col = 0; col < m; ++col) {
                ata[row][col] += v[m - 1 - row] * v[m - 1 - col];
            }
            aty[row] += v[m - 1 - row] * y[i];
        }
    }

    std::vector<double> p2 = solve_linear(ata, aty);
    return polyfit_convert(p2, x);
}

double polyval(const std::vector<double> &p, double x) {
    double y = 0.0;
    for (double coeff : p) {
        y = y * x + coeff;
    }
    return y;
}

std::vector<double> conv_same(const std::vector<double> &x, const std::vector<double> &h) {
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

Matrix deriv1(const Matrix &a, const std::vector<double> &fil) {
    Matrix b(a.rows, a.cols, 0.0);
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

double centroid(const std::vector<double> &x) {
    double sum = std::accumulate(x.begin(), x.end(), 0.0);
    if (sum == 0.0) {
        return 0.0;
    }
    double loc = 0.0;
    for (size_t i = 0; i < x.size(); ++i) {
        loc += (static_cast<double>(i) + 1.0) * x[i];
    }
    return loc / sum;
}

std::vector<double> cent(const std::vector<double> &a, double center) {
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

Image rotate90(const Image &in) {
    Image out(in.cols, in.rows, in.channels, 0.0);
    for (int ch = 0; ch < in.channels; ++ch) {
        for (int r = 0; r < in.rows; ++r) {
            for (int c = 0; c < in.cols; ++c) {
                out.at(in.cols - 1 - c, r, ch) = in.at(r, c, ch);
            }
        }
    }
    return out;
}

struct RotateResult {
    Image image;
    int nlin = 0;
    int npix = 0;
    int rflag = 0;
};

RotateResult rotatev2(const Image &input) {
    RotateResult result{input, input.rows, input.cols, 0};
    int nlin = input.rows;
    int npix = input.cols;
    int mm = (input.channels == 3 || input.channels == 4) ? 1 : 0;

    int nn = 3;
    int row_top = std::max(0, nn - 1);
    int row_bot = std::max(0, nlin - nn - 1);
    int col_left = std::max(0, nn - 1);
    int col_right = std::max(0, npix - nn - 1);

    double mean_bot = 0.0;
    double mean_top = 0.0;
    double mean_right = 0.0;
    double mean_left = 0.0;

    for (int c = 0; c < npix; ++c) {
        mean_bot += input.at(row_bot, c, mm);
        mean_top += input.at(row_top, c, mm);
    }
    mean_bot /= static_cast<double>(npix);
    mean_top /= static_cast<double>(npix);

    for (int r = 0; r < nlin; ++r) {
        mean_right += input.at(r, col_right, mm);
        mean_left += input.at(r, col_left, mm);
    }
    mean_right /= static_cast<double>(nlin);
    mean_left /= static_cast<double>(nlin);

    double testv = std::abs(mean_bot - mean_top);
    double testh = std::abs(mean_right - mean_left);

    if (testv > testh) {
        result.rflag = 1;
        result.image = rotate90(input);
        result.nlin = result.image.rows;
        result.npix = result.image.cols;
    }
    return result;
}

std::vector<double> ahamming(int n, double mid) {
    std::vector<double> data(n, 0.0);
    mid += 0.5;
    double wid1 = mid - 1.0;
    double wid2 = static_cast<double>(n) - mid;
    double wid = std::max(wid1, wid2);
    for (int i = 0; i < n; ++i) {
        double arg = (static_cast<double>(i) + 1.0) - mid;
        data[i] = std::cos(M_PI * arg / wid);
    }
    for (double &v : data) {
        v = 0.54 + 0.46 * v;
    }
    return data;
}

std::vector<double> tukey(int n, double alpha) {
    if (n == 1) {
        return {1.0};
    }
    if (alpha == 0.0) {
        return std::vector<double>(n, 1.0);
    }
    double m = (n - 1) / 2.0;
    std::vector<double> w(n, 0.0);
    for (int k = 0; k <= static_cast<int>(m); ++k) {
        if (k <= alpha * m) {
            w[k] = 0.5 * (1 + std::cos(M_PI * (k / (alpha * m) - 1)));
        } else {
            w[k] = 1.0;
        }
        w[n - 1 - k] = w[k];
    }
    return w;
}

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
        w.resize(n);
        return w;
    }
    std::vector<double> out(n, 0.0);
    int start = static_cast<int>(w.size()) - n;
    for (int i = 0; i < n; ++i) {
        out[i] = w[start + i];
    }
    return out;
}

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

std::vector<double> findedge2(const std::vector<double> &cent, int nlin, int nn) {
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

ProjectResult project2(const Matrix &bb, const std::vector<double> &fitme, int fac) {
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

    int zero_counts = 0;
    for (int i = start; i < start + nn; ++i) {
        if (counts[i - 1] == 0.0) {
            zero_counts++;
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

std::vector<std::complex<double>> fft(const std::vector<std::complex<double>> &x) {
    int n = static_cast<int>(x.size());
    if (n == 1) {
        return x;
    }
    if (n % 2 != 0) {
        std::vector<std::complex<double>> out(n);
        for (int k = 0; k < n; ++k) {
            std::complex<double> sum(0.0, 0.0);
            for (int t = 0; t < n; ++t) {
                double angle = -2.0 * M_PI * k * t / n;
                sum += x[t] * std::complex<double>(std::cos(angle), std::sin(angle));
            }
            out[k] = sum;
        }
        return out;
    }
    std::vector<std::complex<double>> even(n / 2);
    std::vector<std::complex<double>> odd(n / 2);
    for (int i = 0; i < n / 2; ++i) {
        even[i] = x[2 * i];
        odd[i] = x[2 * i + 1];
    }
    even = fft(even);
    odd = fft(odd);
    std::vector<std::complex<double>> out(n);
    for (int k = 0; k < n / 2; ++k) {
        double angle = -2.0 * M_PI * k / n;
        std::complex<double> twiddle(std::cos(angle), std::sin(angle));
        out[k] = even[k] + twiddle * odd[k];
        out[k + n / 2] = even[k] - twiddle * odd[k];
    }
    return out;
}

std::vector<double> findfreq(const Matrix &dat, double val, int imax, int fflag) {
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

Matrix sampeff(const Matrix &dat, const std::vector<double> &val, double del, int fflag) {
    if (dat.rows == 0 || dat.cols < 2) {
        return Matrix();
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
        return Matrix(val.size(), dat.cols - 1, 0.0);
    }

    int nc = dat.cols - 1;
    Matrix eff(static_cast<int>(val.size()), nc, 0.0);
    for (size_t v = 0; v < val.size(); ++v) {
        std::vector<double> freq_sfr = findfreq(dat, val[v], imax, fflag);
        for (int c = 0; c < nc; ++c) {
            double freq = clip_value(freq_sfr[c], 0.0, hs);
            eff(static_cast<int>(v), c) = std::min(std::round(100.0 * freq / hs), 100.0);
        }
    }
    return eff;
}

void rsquare(const std::vector<double> &y, const std::vector<double> &f, double &r2, double &rmse) {
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
    double mean_y = mean(yy);
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

}  // namespace

SfrResult compute_sfr(const Image &input,
                      double del,
                      int npol,
                      int wflag,
                      const std::array<double, 3> &weight) {
    if (input.rows == 0 || input.cols == 0) {
        throw std::runtime_error("Empty input image");
    }
    int nbin = 4;
    double alpha = 1.0;
    npol = std::min(npol, 5);

    Image a = input;
    if (a.channels == 3) {
        Image out(a.rows, a.cols, 4, 0.0);
        for (int r = 0; r < a.rows; ++r) {
            for (int c = 0; c < a.cols; ++c) {
                double lum = weight[0] * a.at(r, c, 0) +
                             weight[1] * a.at(r, c, 1) +
                             weight[2] * a.at(r, c, 2);
                out.at(r, c, 0) = a.at(r, c, 0);
                out.at(r, c, 1) = a.at(r, c, 1);
                out.at(r, c, 2) = a.at(r, c, 2);
                out.at(r, c, 3) = lum;
            }
        }
        a = out;
    }

    if (del > 1.0) {
        del = 25.4 / del;
    }

    RotateResult rot = rotatev2(a);
    a = rot.image;
    int nlin = a.rows;
    int npix = a.cols;
    int ncol = a.channels;

    int left_cols = std::min(5, npix);
    int right_cols = std::min(6, npix);
    double tleft = 0.0;
    double tright = 0.0;
    for (int r = 0; r < nlin; ++r) {
        for (int c = 0; c < left_cols; ++c) {
            tleft += a.at(r, c, 0);
        }
        for (int c = npix - right_cols; c < npix; ++c) {
            tright += a.at(r, c, 0);
        }
    }

    std::vector<double> fil1 = {0.5, -0.5};
    std::vector<double> fil2 = {0.5, 0.0, -0.5};
    if (tleft > tright) {
        fil1 = {-0.5, 0.5};
        fil2 = {-0.5, 0.0, 0.5};
    }

    std::vector<double> win1;
    if (wflag != 0) {
        win1 = ahamming(npix, (npix + 1) / 2.0);
    } else {
        win1 = tukey2(npix, alpha, (npix + 1) / 2.0);
        for (double &v : win1) {
            v = 0.95 * v + 0.05;
        }
    }

    std::vector<std::vector<double>> loc(ncol, std::vector<double>(nlin, 0.0));
    std::vector<std::vector<double>> fitme(ncol);
    std::vector<std::vector<double>> fitme1(ncol);

    for (int color = 0; color < ncol; ++color) {
        Matrix plane(nlin, npix, 0.0);
        for (int r = 0; r < nlin; ++r) {
            for (int c = 0; c < npix; ++c) {
                plane(r, c) = a.at(r, c, color);
            }
        }
        Matrix deriv = deriv1(plane, fil1);
        for (int n = 0; n < nlin; ++n) {
            std::vector<double> row(npix, 0.0);
            for (int c = 0; c < npix; ++c) {
                row[c] = deriv(n, c) * win1[c];
            }
            loc[color][n] = centroid(row) - 0.5;
        }
        fitme[color] = findedge2(loc[color], nlin, npol);

        for (int n = 0; n < nlin; ++n) {
            double place = polyval(fitme[color], static_cast<double>(n));
            std::vector<double> win2 = (wflag != 0) ? ahamming(npix, place)
                                                    : tukey2(npix, alpha, place);
            if (wflag == 0) {
                for (double &v : win2) {
                    v = 0.95 * v + 0.05;
                }
            }
            std::vector<double> row(npix, 0.0);
            for (int c = 0; c < npix; ++c) {
                row[c] = deriv(n, c) * win2[c];
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
        nlin1 = static_cast<int>(std::round(std::floor(nlin * std::abs(slope_ref)) / std::abs(slope_ref)));
    }
    nlin1 = std::max(1, std::min(nlin1, nlin));
    if (nlin1 < nlin) {
        Image cropped(nlin1, npix, ncol, 0.0);
        for (int r = 0; r < nlin1; ++r) {
            for (int c = 0; c < npix; ++c) {
                for (int ch = 0; ch < ncol; ++ch) {
                    cropped.at(r, c, ch) = a.at(r, c, ch);
                }
            }
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
        for (double &m : misreg) {
            m *= delfac;
        }
    }

    int nn = static_cast<int>(std::ceil(npix * nbin));
    int nn2 = static_cast<int>(std::floor(nn / 2.0)) + 1;
    std::vector<double> dcorr = fir2fix(nn2, 3);
    int freqlim = (nbin == 1) ? 2 : 1;
    int nn2out = static_cast<int>(std::round(nn2 * freqlim / 2.0));

    Matrix mtf(nn, ncol, 0.0);
    std::vector<double> esf_last;

    for (int color = 0; color < ncol; ++color) {
        Matrix plane(nlin, npix, 0.0);
        for (int r = 0; r < nlin; ++r) {
            for (int c = 0; c < npix; ++c) {
                plane(r, c) = a.at(r, c, color);
            }
        }
        ProjectResult proj = project2(plane, fitme[color], nbin);
        esf_last = proj.point;
        Matrix esf_mat(1, static_cast<int>(esf_last.size()), 0.0);
        for (int i = 0; i < static_cast<int>(esf_last.size()); ++i) {
            esf_mat(0, i) = esf_last[i];
        }
        Matrix deriv = deriv1(esf_mat, fil2);
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
        std::vector<double> win = (wflag != 0) ? ahamming(nn, center) : tukey2(nn, alpha, center);
        for (int i = 0; i < nn; ++i) {
            c[i] *= win[i];
        }

        std::vector<std::complex<double>> cx(nn);
        for (int i = 0; i < nn; ++i) {
            cx[i] = std::complex<double>(c[i], 0.0);
        }
        std::vector<std::complex<double>> fx = fft(cx);
        double dc0 = std::abs(fx[0]);
        for (int i = 0; i < nn2; ++i) {
            double val = (dc0 == 0.0) ? 0.0 : std::abs(fx[i]) / dc0;
            val *= dcorr[i];
            mtf(i, color) = val;
        }
    }

    std::vector<double> freq(nn, 0.0);
    for (int n = 0; n < nn; ++n) {
        freq[n] = static_cast<double>(n) / (del2 * nn);
    }
    Matrix dat(nn2out, ncol + 1, 0.0);
    for (int i = 0; i < nn2out; ++i) {
        dat(i, 0) = freq[i];
        for (int c = 0; c < ncol; ++c) {
            dat(i, c + 1) = mtf(i, c);
        }
    }

    int fit_cols = static_cast<int>(fitme[0].size());
    Matrix fitout(ncol, (ncol > 2) ? fit_cols + 1 : fit_cols, 0.0);
    for (int r = 0; r < ncol; ++r) {
        for (int c = 0; c < fit_cols; ++c) {
            fitout(r, c) = fitme[r][c];
        }
        if (ncol > 2) {
            fitout(r, fit_cols) = misreg[r];
        }
    }

    std::vector<double> val = {0.1, 0.5};
    Matrix eff = sampeff(dat, val, delimage, 0);
    std::vector<double> freq_sfr = findfreq(dat, 0.5, dat.rows, 0);
    double sfr50 = freq_sfr.empty() ? 0.0 : freq_sfr[0];

    SfrResult result;
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

}  // namespace sfrmat5
