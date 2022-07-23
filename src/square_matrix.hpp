/**
 * @file square_matrix.hpp
 * @author remzerrr (remi.helleboid@gmail.com)
 * @brief
 * @version 0.1
 * @date 2022-07-23
 *
 * @copyright Copyright (c) 2022
 *
 */

#pragma once

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <memory>
#include <numeric>
#include <random>
#include <vector>

class Matrix {
 private:
    std::size_t         m_size;
    std::vector<double> m_data;

 public:
    Matrix(std::size_t size) : m_size(size), m_data(size * size) {}
    Matrix(std::size_t size, double value) : m_size(size), m_data(size * size, value) {}
    Matrix(std::size_t size, std::vector<double> data) : m_size(size), m_data(data) {}
    Matrix(const Matrix&) = default;
    Matrix(Matrix&&)      = default;
    Matrix& operator=(const Matrix&) = default;
    Matrix& operator=(Matrix&&) = default;
    ~Matrix()                   = default;

    std::size_t size() const { return m_size; }
    double&     operator()(std::size_t i, std::size_t j) { return m_data[i * m_size + j]; }
    double      operator()(std::size_t i, std::size_t j) const { return m_data[i * m_size + j]; }

    void set_identity() {
        std::fill(m_data.begin(), m_data.end(), 0);
        for (std::size_t i = 0; i < m_size; ++i) {
            (*this)(i, i) = 1;
        }
    }

    void set_zero() { std::fill(m_data.begin(), m_data.end(), 0); }

    void set_value(double value) { std::fill(m_data.begin(), m_data.end(), value); }

    void set_random(const std::string& distribution, double param1, double param2) {
        std::random_device rd;
        std::mt19937       gen(rd());
        if (distribution == "uniform") {
            std::uniform_real_distribution<> dis(param1, param2);
            std::generate(m_data.begin(), m_data.end(), [&]() { return dis(gen); });
        } else if (distribution == "normal") {
            std::normal_distribution<> dis(param1, param2);
            std::generate(m_data.begin(), m_data.end(), [&]() { return dis(gen); });
        } else if (distribution == "exponential") {
            std::exponential_distribution<> dis(param1);
            std::generate(m_data.begin(), m_data.end(), [&]() { return dis(gen); });
        } else if (distribution == "poisson") {
            std::poisson_distribution<> dis(param1);
            std::generate(m_data.begin(), m_data.end(), [&]() { return dis(gen); });
        } else {
            throw std::runtime_error("Unknown distribution");
        }
    }

    Matrix& operator+=(const Matrix& rhs) {
        assert(m_size == rhs.m_size);
        std::transform(m_data.begin(), m_data.end(), rhs.m_data.begin(), m_data.begin(), std::plus<double>());
        return *this;
    }

    Matrix& operator-=(const Matrix& rhs) {
        assert(m_size == rhs.m_size);
        std::transform(m_data.begin(), m_data.end(), rhs.m_data.begin(), m_data.begin(), std::minus<double>());
        return *this;
    }

    Matrix& operator*=(const Matrix& rhs) {
        assert(m_size == rhs.m_size);
        std::vector<double> tmp(m_size * m_size);
        for (std::size_t i = 0; i < m_size; ++i) {
            for (std::size_t j = 0; j < m_size; ++j) {
                tmp[i * m_size + j] = std::inner_product(m_data.begin() + i * m_size,
                                                         m_data.begin() + (i + 1) * m_size,
                                                         rhs.m_data.begin() + j * m_size,
                                                         0.0);
            }
        }
        m_data = tmp;
        return *this;
    }

    Matrix& operator*=(double rhs) {
        std::transform(m_data.begin(), m_data.end(), m_data.begin(), [rhs](double x) { return x * rhs; });
        return *this;
    }

    Matrix& operator/=(double rhs) {
        std::transform(m_data.begin(), m_data.end(), m_data.begin(), [rhs](double x) { return x / rhs; });
        return *this;
    }

    Matrix operator+(const Matrix& rhs) const {
        assert(m_size == rhs.m_size);
        Matrix tmp(*this);
        tmp += rhs;
        return tmp;
    }

    Matrix operator-(const Matrix& rhs) const {
        assert(m_size == rhs.m_size);
        Matrix tmp(*this);
        tmp -= rhs;
        return tmp;
    }

    Matrix operator*(const Matrix& rhs) const {
        assert(m_size == rhs.m_size);
        Matrix tmp(*this);
        tmp *= rhs;
        return tmp;
    }

    Matrix operator*(double rhs) const {
        Matrix tmp(*this);
        tmp *= rhs;
        return tmp;
    }

    Matrix operator/(double rhs) const {
        Matrix tmp(*this);
        tmp /= rhs;
        return tmp;
    }

    Matrix operator-() const {
        Matrix tmp(*this);
        std::transform(tmp.m_data.begin(), tmp.m_data.end(), tmp.m_data.begin(), [](double x) { return -x; });
        return tmp;
    }

    bool operator==(const Matrix& rhs) const {
        assert(m_size == rhs.m_size);
        return std::equal(m_data.begin(), m_data.end(), rhs.m_data.begin());
    }

    bool operator!=(const Matrix& rhs) const {
        assert(m_size == rhs.m_size);
        return !std::equal(m_data.begin(), m_data.end(), rhs.m_data.begin());
    }

    Matrix exponential() const {
        Matrix tmp(*this);
        std::transform(tmp.m_data.begin(), tmp.m_data.end(), tmp.m_data.begin(), [](double x) { return std::exp(x); });
        return tmp;
    }

    double norm() const { return std::sqrt(std::inner_product(m_data.begin(), m_data.end(), m_data.begin(), 0.0)); }

    double norm_squared() const { return std::inner_product(m_data.begin(), m_data.end(), m_data.begin(), 0.0); }

    double min() const { return *std::min_element(m_data.begin(), m_data.end()); }

    double max() const { return *std::max_element(m_data.begin(), m_data.end()); }

    double sum() const { return std::accumulate(m_data.begin(), m_data.end(), 0.0); }

    bool is_diagonal_dominant() const {
        for (std::size_t i = 0; i < m_size; ++i) {
            double sum = 0;
            for (std::size_t j = 0; j < m_size; ++j) {
                if (i != j) {
                    sum += std::abs(m_data[i * m_size + j]);
                }
            }
            if (sum > m_data[i * m_size + i]) {
                return false;
            }
        }
        return true;
    }

    bool is_symmetric() const {
        for (std::size_t i = 0; i < m_size; ++i) {
            for (std::size_t j = 0; j < m_size; ++j) {
                if (m_data[i * m_size + j] != m_data[j * m_size + i]) {
                    return false;
                }
            }
        }
        return true;
    }

    bool is_diagonal() const {
        for (std::size_t i = 0; i < m_size; ++i) {
            for (std::size_t j = 0; j < m_size; ++j) {
                if (i != j && m_data[i * m_size + j] != 0) {
                    return false;
                }
            }
        }
        return true;
    }

    Matrix transpose() const {
        Matrix tmp(m_size);
        for (std::size_t i = 0; i < m_size; ++i) {
            for (std::size_t j = 0; j < m_size; ++j) {
                tmp(i, j) = m_data[j * m_size + i];
            }
        }
        return tmp;
    }

    Matrix& transpose_in_place() {
        for (std::size_t i = 0; i < m_size; ++i) {
            for (std::size_t j = 0; j < m_size; ++j) {
                std::swap(m_data[i * m_size + j], m_data[j * m_size + i]);
            }
        }
        return *this;
    }

    double trace() const {
        double trace = 0;
        for (std::size_t i = 0; i < m_size; ++i) {
            trace += m_data[i * m_size + i];
        }
        return trace;
    }

    Matrix sub_matrix_rank_1(std::size_t row_to_removed, std::size_t col_to_removed) const {
        assert(row_to_removed < m_size && col_to_removed < m_size);
        Matrix      tmp(m_size - 1);
        std::size_t row_idx = 0;
        for (std::size_t i = 0; i < m_size; ++i) {
            if (i != row_to_removed) {
                std::size_t col_idx = 0;
                for (std::size_t j = 0; j < m_size; ++j) {
                    if (j != col_to_removed) {
                        tmp(row_idx, col_idx) = m_data[i * m_size + j];
                        ++col_idx;
                    }
                }
                ++row_idx;
            }
        }
        return tmp;
    }

    double compute_determinant() const {
        double det = 0.0;
        if (m_size == 1) {
            det = m_data[0];
        } else if (m_size == 2) {
            det = m_data[0] * m_data[3] - m_data[1] * m_data[2];
        } else {
            for (std::size_t i = 0; i < m_size; ++i) {
                det += m_data[i * m_size + i] * sub_matrix_rank_1(i, 0).compute_determinant();
            }
        }
        return det;
    }

    char* serialize() const {
        char* buffer = new char[m_size * m_size * sizeof(double)];
        std::memcpy(buffer, m_data.data(), m_size * m_size * sizeof(double));
        return buffer;
    }

    friend Matrix deserialize(const char* buffer) {
        std::size_t size = sqrt(std::strlen(buffer) / sizeof(double));
        Matrix      tmp(size);
        std::memcpy(tmp.m_data.data(), buffer, size * size * sizeof(double));
        return tmp;
    }

    void export_bitmap(const std::string& filename) const {
        std::ofstream file(filename, std::ios::out | std::ios::binary);
        if (!file) {
            throw std::runtime_error("Could not open file " + filename);
        }
        file << "P6\n" << m_size << " " << m_size << "\n255\n";
        for (std::size_t i = 0; i < m_size; ++i) {
            for (std::size_t j = 0; j < m_size; ++j) {
                file << (unsigned char)(255 * m_data[i * m_size + j]);
            }
        }
        file.close();
    }

    friend std::ostream& operator<<(std::ostream& os, const Matrix& m) {
        for (std::size_t i = 0; i < m.m_size; ++i) {
            for (std::size_t j = 0; j < m.m_size; ++j) {
                os << std::setw(10) << m(i, j) << " ";
            }
            os << std::endl;
        }
        return os;
    }

    void print() const {
        for (std::size_t i = 0; i < m_size; ++i) {
            for (std::size_t j = 0; j < m_size; ++j) {
                std::cout << std::setw(16) << m_data[i * m_size + j];
            }
            std::cout << std::endl;
        }
    }

    void export_to_file(const std::string& filename) const {
        std::ofstream ofs(filename);
        ofs << *this;
        ofs.close();
    }
};