#include "sfrmat5.h"

#include <cmath>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"

namespace {

using Scalar = double;

struct Image {
    int rows = 0;
    int cols = 0;
    int channels = 0;
    std::vector<sfrmat5::Matrix<Scalar>> planes;

    /// Constructs an empty test image container.
    Image() = default;

    /// Constructs a test image with channel planes initialized to a constant value.
    Image(int r, int c, int ch, Scalar value = static_cast<Scalar>(0))
        : rows(r), cols(c), channels(ch), planes(ch, sfrmat5::Matrix<Scalar>(r, c)) {
        for (int i = 0; i < ch; ++i) {
            planes[i].setConstant(value);
        }
    }
};

/// Returns true when a value is numerically close to zero.
bool nearly_zero(double v) {
    return std::abs(v) < 1e-9;
}

/// Returns true when two finite values are within the requested tolerance.
bool nearly_equal(double actual, double expected, double tol) {
    return std::isfinite(actual) && std::abs(actual - expected) <= tol;
}

/// Verifies that the first SFR data column is a strictly increasing frequency axis.
bool check_frequency_axis(const sfrmat5::Matrix<Scalar>& dat) {
    if (dat.rows() < 2 || dat.cols() < 2) {
        return false;
    }
    double prev = dat(0, 0);
    for (int i = 1; i < dat.rows(); ++i) {
        double cur = dat(i, 0);
        if (!(cur > prev)) {
            return false;
        }
        prev = cur;
    }
    return true;
}

/// Checks a scalar value against an expected value and reports failures.
bool check_value(const char* label, double actual, double expected, double tol) {
    if (nearly_equal(actual, expected, tol)) {
        return true;
    }
    std::cerr << label << " mismatch: expected " << expected << ", got " << actual
              << ", tolerance " << tol << "\n";
    return false;
}

/// Checks one matrix element against an expected value and reports failures.
bool check_matrix_value(const char* label, const sfrmat5::Matrix<Scalar>& m, int row, int col,
                        double expected, double tol) {
    if (row >= m.rows() || col >= m.cols()) {
        std::cerr << label << " index out of range at (" << row << ", " << col << ")\n";
        return false;
    }
    return check_value(label, m(row, col), expected, tol);
}

/// Loads an image file into planar scalar channels using stb_image.
Image load_image(const std::string& path) {
    int width = 0;
    int height = 0;
    int source_channels = 0;
    if (!stbi_info(path.c_str(), &width, &height, &source_channels)) {
        throw std::runtime_error(std::string("Failed to inspect image: ") + stbi_failure_reason());
    }

    const int output_channels = (source_channels == 1) ? 1 : 3;
    std::unique_ptr<unsigned char, decltype(&stbi_image_free)> data(
        stbi_load(path.c_str(), &width, &height, &source_channels, output_channels),
        stbi_image_free);
    if (!data) {
        throw std::runtime_error(std::string("Failed to load image: ") + stbi_failure_reason());
    }

    Image img(height, width, output_channels, static_cast<Scalar>(0));
    for (int row = 0; row < img.rows; ++row) {
        for (int col = 0; col < img.cols; ++col) {
            const size_t pixel_offset =
                (static_cast<size_t>(row) * img.cols + static_cast<size_t>(col)) * img.channels;
            for (int ch = 0; ch < img.channels; ++ch) {
                img.planes[ch](row, col) = static_cast<Scalar>(data.get()[pixel_offset + ch]);
            }
        }
    }
    return img;
}

/// Flattens image planes into the public SfrMat5 planar pixel layout.
std::vector<Scalar> extract_planar_pixels(const Image& img) {
    std::vector<Scalar> pixels(static_cast<size_t>(img.rows) * static_cast<size_t>(img.cols) *
                               static_cast<size_t>(img.channels));
    const size_t plane_size = static_cast<size_t>(img.rows) * static_cast<size_t>(img.cols);
    for (int ch = 0; ch < img.channels; ++ch) {
        const size_t channel_offset = static_cast<size_t>(ch) * plane_size;
        for (int row = 0; row < img.rows; ++row) {
            const size_t row_offset = channel_offset + static_cast<size_t>(row) * img.cols;
            for (int col = 0; col < img.cols; ++col) {
                pixels[row_offset + col] = img.planes[ch](row, col);
            }
        }
    }
    return pixels;
}

} // namespace

/// Runs the regression test against the example edge image.
int main() {
    std::string path = "Example_Images/Test_edge1.bmp";
    Image img = load_image(path);
    auto pixels = std::make_unique<std::vector<Scalar>>(extract_planar_pixels(img));
    sfrmat5::SfrMat5<Scalar> sfr;
    sfrmat5::SfrResult<Scalar> result =
        sfr.compute(std::move(pixels), img.cols, img.rows, img.channels);

    if (result.dat.rows() == 0 || result.dat.cols() < 2) {
        std::cerr << "SFR data missing\n";
        return 1;
    }
    if (result.dat.rows() != 125 || result.dat.cols() != 5) {
        std::cerr << "Unexpected SFR data dimensions: " << result.dat.rows() << "x"
                  << result.dat.cols() << "\n";
        return 1;
    }
    if (!check_frequency_axis(result.dat)) {
        std::cerr << "Frequency axis not increasing\n";
        return 1;
    }
    if (nearly_zero(result.sfr50) || std::isnan(result.sfr50)) {
        std::cerr << "SFR50 invalid\n";
        return 1;
    }
    if (result.e.rows() == 0 || result.e.cols() == 0) {
        std::cerr << "Sampling efficiency missing\n";
        return 1;
    }
    if (result.e.rows() != 2 || result.e.cols() != 4) {
        std::cerr << "Unexpected sampling efficiency dimensions: " << result.e.rows() << "x"
                  << result.e.cols() << "\n";
        return 1;
    }

    const double value_tol = 1e-5;
    const double freq_tol = 1e-6;
    bool numerical_ok = true;
    numerical_ok &= check_value("SFR50", result.sfr50, 0.269805, value_tol);
    numerical_ok &= check_value("del2", result.del2, 0.248855, value_tol);

    numerical_ok &= check_matrix_value("sampling efficiency 10% R", result.e, 0, 0, 85.0, 0.0);
    numerical_ok &= check_matrix_value("sampling efficiency 10% G", result.e, 0, 1, 85.0, 0.0);
    numerical_ok &= check_matrix_value("sampling efficiency 10% B", result.e, 0, 2, 86.0, 0.0);
    numerical_ok &= check_matrix_value("sampling efficiency 10% L", result.e, 0, 3, 85.0, 0.0);
    numerical_ok &= check_matrix_value("sampling efficiency 50% R", result.e, 1, 0, 55.0, 0.0);
    numerical_ok &= check_matrix_value("sampling efficiency 50% G", result.e, 1, 1, 55.0, 0.0);
    numerical_ok &= check_matrix_value("sampling efficiency 50% B", result.e, 1, 2, 56.0, 0.0);
    numerical_ok &= check_matrix_value("sampling efficiency 50% L", result.e, 1, 3, 55.0, 0.0);

    numerical_ok &= check_matrix_value("dat[0,0]", result.dat, 0, 0, 0.0, freq_tol);
    numerical_ok &= check_matrix_value("dat[0,1]", result.dat, 0, 1, 1.0, value_tol);
    numerical_ok &= check_matrix_value("dat[1,0]", result.dat, 1, 0, 0.00810162, freq_tol);
    numerical_ok &= check_matrix_value("dat[1,4]", result.dat, 1, 4, 0.994307, value_tol);
    numerical_ok &= check_matrix_value("dat[33,0]", result.dat, 33, 0, 0.267353, freq_tol);
    numerical_ok &= check_matrix_value("dat[33,4]", result.dat, 33, 4, 0.511463, value_tol);
    numerical_ok &= check_matrix_value("dat[61,0]", result.dat, 61, 0, 0.494199, freq_tol);
    numerical_ok &= check_matrix_value("dat[61,4]", result.dat, 61, 4, 0.0146672, value_tol);
    numerical_ok &= check_matrix_value("dat[124,0]", result.dat, 124, 0, 1.0046, value_tol);
    numerical_ok &= check_matrix_value("dat[124,4]", result.dat, 124, 4, 0.0797956, value_tol);
    if (!numerical_ok) {
        return 1;
    }

    std::cout << "sfrmat5 basic test passed\n";
    std::cout << "SFR50: " << result.sfr50 << "\n";
    if (result.e.rows() > 0 && result.e.cols() > 0) {
        std::cout << "Sampling efficiency (10%): ";
        for (int c = 0; c < result.e.cols(); ++c) {
            std::cout << result.e(0, c);
            if (c + 1 < result.e.cols()) {
                std::cout << ", ";
            }
        }
        std::cout << "\n";
    }
    std::cout << "First " << result.dat.rows() << " SFR rows (freq, mtf...):\n";
    for (int i = 0; i < result.dat.rows(); ++i) {
        for (int c = 0; c < result.dat.cols(); ++c) {
            std::cout << result.dat(i, c);
            if (c + 1 < result.dat.cols()) {
                std::cout << ", ";
            }
        }
        std::cout << "\n";
    }
    return 0;
}
