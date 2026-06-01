#include "sfrmat5.h"

#include <cmath>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace {

using Scalar = double;

struct BmpImage {
    int rows = 0;
    int cols = 0;
    int channels = 0;
    std::vector<sfrmat5::Matrix<Scalar>> planes;

    BmpImage() = default;
    BmpImage(int r, int c, int ch, Scalar value = static_cast<Scalar>(0))
        : rows(r), cols(c), channels(ch), planes(ch, sfrmat5::Matrix<Scalar>(r, c)) {
        for (int i = 0; i < ch; ++i) {
            planes[i].setConstant(value);
        }
    }
};

bool nearly_zero(double v) {
    return std::abs(v) < 1e-9;
}

bool nearly_equal(double actual, double expected, double tol) {
    return std::isfinite(actual) && std::abs(actual - expected) <= tol;
}

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

bool check_value(const char* label, double actual, double expected, double tol) {
    if (nearly_equal(actual, expected, tol)) {
        return true;
    }
    std::cerr << label << " mismatch: expected " << expected << ", got " << actual
              << ", tolerance " << tol << "\n";
    return false;
}

bool check_matrix_value(const char* label, const sfrmat5::Matrix<Scalar>& m, int row, int col,
                        double expected, double tol) {
    if (row >= m.rows() || col >= m.cols()) {
        std::cerr << label << " index out of range at (" << row << ", " << col << ")\n";
        return false;
    }
    return check_value(label, m(row, col), expected, tol);
}

uint16_t read_u16(std::ifstream& in) {
    uint8_t b0 = 0;
    uint8_t b1 = 0;
    in.read(reinterpret_cast<char*>(&b0), 1);
    in.read(reinterpret_cast<char*>(&b1), 1);
    return static_cast<uint16_t>(b0 | (b1 << 8));
}

uint32_t read_u32(std::ifstream& in) {
    uint8_t b[4] = {0, 0, 0, 0};
    in.read(reinterpret_cast<char*>(b), 4);
    return static_cast<uint32_t>(b[0] | (b[1] << 8) | (b[2] << 16) | (b[3] << 24));
}

int32_t read_i32(std::ifstream& in) {
    return static_cast<int32_t>(read_u32(in));
}

BmpImage load_bmp(const std::string& path) {
    std::ifstream in(path, std::ios::binary);
    if (!in) {
        throw std::runtime_error("Failed to open BMP");
    }

    uint16_t bfType = read_u16(in);
    if (bfType != 0x4D42) {
        throw std::runtime_error("Not a BMP file");
    }
    uint32_t bfSize = read_u32(in);
    (void)bfSize;
    read_u16(in);
    read_u16(in);
    uint32_t bfOffBits = read_u32(in);

    uint32_t biSize = read_u32(in);
    if (biSize < 40) {
        throw std::runtime_error("Unsupported BMP header");
    }
    int32_t width = read_i32(in);
    int32_t height = read_i32(in);
    uint16_t planes = read_u16(in);
    uint16_t bitCount = read_u16(in);
    uint32_t compression = read_u32(in);
    read_u32(in);
    read_i32(in);
    read_i32(in);
    read_u32(in);
    read_u32(in);

    if (planes != 1 || (bitCount != 8 && bitCount != 24) || compression != 0) {
        throw std::runtime_error("Unsupported BMP format");
    }

    if (biSize > 40) {
        in.seekg(static_cast<std::streamoff>(biSize - 40), std::ios::cur);
    }

    bool bottom_up = true;
    if (height < 0) {
        bottom_up = false;
        height = -height;
    }

    if (bitCount == 8) {
        int palette_entries = static_cast<int>((bfOffBits - 54) / 4);
        in.seekg(54 + palette_entries * 4, std::ios::beg);
    } else {
        in.seekg(static_cast<std::streamoff>(bfOffBits), std::ios::beg);
    }

    int rows = height;
    int cols = width;
    int channels = (bitCount == 24) ? 3 : 1;
    BmpImage img(rows, cols, channels, static_cast<Scalar>(0));

    int row_bytes = ((bitCount * cols + 31) / 32) * 4;
    std::vector<uint8_t> row(row_bytes, 0);
    for (int r = 0; r < rows; ++r) {
        int dst_row = bottom_up ? (rows - 1 - r) : r;
        in.read(reinterpret_cast<char*>(row.data()), row_bytes);
        if (!in) {
            throw std::runtime_error("BMP read failed");
        }
        if (bitCount == 24) {
            for (int c = 0; c < cols; ++c) {
                int idx = c * 3;
                uint8_t b = row[idx];
                uint8_t g = row[idx + 1];
                uint8_t rch = row[idx + 2];
                img.planes[0](dst_row, c) = static_cast<Scalar>(rch);
                img.planes[1](dst_row, c) = static_cast<Scalar>(g);
                img.planes[2](dst_row, c) = static_cast<Scalar>(b);
            }
        } else {
            for (int c = 0; c < cols; ++c) {
                img.planes[0](dst_row, c) = static_cast<Scalar>(row[c]);
            }
        }
    }
    return img;
}

std::vector<Scalar> extract_planar_pixels(const BmpImage& img) {
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

int main() {
    std::string path = "Example_Images/Test_edge1.bmp";
    BmpImage img = load_bmp(path);
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
