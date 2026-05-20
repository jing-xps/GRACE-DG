#include "MathUtils.hpp"
#include "Globals.hpp"
#include <cmath>
#include <algorithm>

// This file contains utility functions for mathematical operations
// including Legendre polynomials and their derivatives
// as well as initialization of quadrature points and weights.

/**
 * @brief Computes Legendre polynomial P_i(x) for i <= 4 (unrolled version)
 * @param i Polynomial order (0–4)
 * @param x Evaluation point
 * @return P_i(x)
 */
double legendre(int i, double x) noexcept {
    switch (i) {
        case 0: return 1.0;
        case 1: return x;
        case 2: return 0.5 * (3.0 * x * x - 1.0);
        case 3: return 0.5 * (5.0 * x * x * x - 3.0 * x);
        case 4: return (35.0 * x * x * x * x - 30.0 * x * x + 3.0) / 8.0;
        default: return 0.0;
    }
}

/**
 * @brief Computes derivative of Legendre polynomial P_i'(x) for i <= 4 (unrolled version)
 * @param i Polynomial order (0–4)
 * @param x Evaluation point
 * @return P_i'(x)
 */
double Dlegendre(int i, double x) noexcept {
    switch (i) {
        case 0: return 0.0;
        case 1: return 1.0;
        case 2: return 3.0 * x;
        case 3: return 0.5 * (15.0 * x * x - 3.0);
        case 4: return 0.5 * (35.0 * x * x * x - 15.0 * x);
        default: return 0.0;
    }
}

/**
 * @brief Initialize precomputed Legendre polynomial table laqp[alpha][i]
 * @return 2D vector of size quad_size x (POLY_ORDER + 1)
 */
std::vector<std::vector<double>> init_laqp() {
    std::vector<std::vector<double>> laqp(Q_SIZE, std::vector<double>(POLY_ORDER + 1));
    for (int alpha = 0; alpha < Q_SIZE; ++alpha)
        for (int i = 0; i <= POLY_ORDER; ++i)
            laqp[alpha][i] = legendre(i, Q_POINTS[alpha]);
    return laqp;
}

/**
 *@brief Initialize precomputed Legendre polynomial table Dlaqp[alpha][i]
 * @return 2D vector of size quad_size x (POLY_ORDER + 1)
 */
std::vector<std::vector<double>> init_Dlaqp() {
    std::vector<std::vector<double>> Dlaqp(Q_SIZE, std::vector<double>(POLY_ORDER + 1));
    for (int alpha = 0; alpha < Q_SIZE; ++alpha)
        for (int i = 0; i <= POLY_ORDER; ++i)
            Dlaqp[alpha][i] = Dlegendre(i, Q_POINTS[alpha]);
    return Dlaqp;
}

/**
 * @brief Finds the bin index for a given Value
 * @param cell Bin boundaries
 * @param val Value to locate
 * @return Bin index
 */
int find_bin_fast(const std::vector<double>& cell, double val) noexcept {
    auto it = std::lower_bound(cell.begin(), cell.end(), val);
    return static_cast<int>(it - cell.begin()) - 1;
}
