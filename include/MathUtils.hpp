#pragma once
#include <vector>

/// @brief Computes Legendre polynomial P_i(x) for i <= 4
double legendre(int i, double x) noexcept;

/// @brief Computes derivative of Legendre polynomial P_i'(x) for i <= 4
double Dlegendre(int i, double x) noexcept;

/// @brief Initialize precomputed Legendre polynomial table laqp
std::vector<std::vector<double>> init_laqp();

/// @brief Initialize precomputed Legendre polynomial table Dlaqp
std::vector<std::vector<double>> init_Dlaqp();

/// @brief Finds the bin index for a given value fast
int find_bin_fast(const std::vector<double>& cell, double val) noexcept;
