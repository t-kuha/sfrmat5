#include "sfrmat5.h"

#include <cstdint>
#include <cmath>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>

namespace {

using Scalar = double;

bool nearly_zero(double v) {
    return std::abs(v) < 1e-9;
}

bool check_frequency_axis(const sfrmat5::Matrix<Scalar> &dat) {
    if (dat.rows < 2 || dat.cols < 2) {
        return false;
    }
    double prev = dat(0, 0);
    for (int i = 1; i < dat.rows; ++i) {
        double cur = dat(i, 0);
        if (!(cur > prev)) {
            return false;
        }
        prev = cur;
    }
    return true;
}

uint16_t read_u16(std::ifstream &in) {
    uint8_t b0 = 0;
    uint8_t b1 = 0;
    in.read(reinterpret_cast<char *>(&b0), 1);
    in.read(reinterpret_cast<char *>(&b1), 1);
    return static_cast<uint16_t>(b0 | (b1 << 8));
}

uint32_t read_u32(std::ifstream &in) {
    uint8_t b[4] = {0, 0, 0, 0};
    in.read(reinterpret_cast<char *>(b), 4);
    return static_cast<uint32_t>(b[0] | (b[1] << 8) | (b[2] << 16) | (b[3] << 24));
}

int32_t read_i32(std::ifstream &in) {
    return static_cast<int32_t>(read_u32(in));
}

sfrmat5::Image<Scalar> load_bmp(const std::string &path) {
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
    sfrmat5::Image<Scalar> img(rows, cols, channels, static_cast<Scalar>(0));

    int row_bytes = ((bitCount * cols + 31) / 32) * 4;
    std::vector<uint8_t> row(row_bytes, 0);
    for (int r = 0; r < rows; ++r) {
        int dst_row = bottom_up ? (rows - 1 - r) : r;
        in.read(reinterpret_cast<char *>(row.data()), row_bytes);
        if (!in) {
            throw std::runtime_error("BMP read failed");
        }
        if (bitCount == 24) {
            for (int c = 0; c < cols; ++c) {
                int idx = c * 3;
                uint8_t b = row[idx];
                uint8_t g = row[idx + 1];
                uint8_t rch = row[idx + 2];
                img.at(dst_row, c, 0) = static_cast<Scalar>(rch);
                img.at(dst_row, c, 1) = static_cast<Scalar>(g);
                img.at(dst_row, c, 2) = static_cast<Scalar>(b);
            }
        } else {
            for (int c = 0; c < cols; ++c) {
                img.at(dst_row, c, 0) = static_cast<Scalar>(row[c]);
            }
        }
    }
    return img;
}

}  // namespace

int main() {
    std::string path = "Example_Images/Test_edge1.bmp";
    sfrmat5::Image<Scalar> img = load_bmp(path);
    sfrmat5::SfrMat5<Scalar> sfr;
    sfrmat5::SfrResult<Scalar> result = sfr.compute(img);

    if (result.dat.rows == 0 || result.dat.cols < 2) {
        std::cerr << "SFR data missing\n";
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
    if (result.e.rows == 0 || result.e.cols == 0) {
        std::cerr << "Sampling efficiency missing\n";
        return 1;
    }

    std::cout << "sfrmat5 basic test passed\n";
    std::cout << "SFR50: " << result.sfr50 << "\n";
    if (result.e.rows > 0 && result.e.cols > 0) {
        std::cout << "Sampling efficiency (10%): ";
        for (int c = 0; c < result.e.cols; ++c) {
            std::cout << result.e(0, c);
            if (c + 1 < result.e.cols) {
                std::cout << ", ";
            }
        }
        std::cout << "\n";
    }
    std::cout << "First " << result.dat.rows << " SFR rows (freq, mtf...):\n";
    for (int i = 0; i < result.dat.rows; ++i) {
        for (int c = 0; c < result.dat.cols; ++c) {
            std::cout << result.dat(i, c);
            if (c + 1 < result.dat.cols) {
                std::cout << ", ";
            }
        }
        std::cout << "\n";
    }
    return 0;
}
