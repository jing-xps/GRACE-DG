#include "DGSolver.hpp"
#include "Globals.hpp"
#include "MathUtils.hpp"
#include "Physics.hpp"
#include <cmath>
#include <iostream>
#include <limits>
#include <algorithm>

// This file implements the core functions for the Discontinuous Galerkin solver
// including initialization of mass bins and density distribution, calculation of cell averages, slope limiters, reconstruction
// and flux calculations for coagulation and fragmentation processes.

// ==========================================================================
// Initialization functions for mass bins and density distribution
// =========================================================================

/**
 * @brief Initialize mass bins with logarithmic spacing
 * @param h Output: Cell widths
 * @param x Output: Cell centers
 * @param cell Output: Cell boundaries
 */
void initialize_mass_bin(std::vector<double>& h, std::vector<double>& x, std::vector<double>& cell) {
    double dx = (std::log10(X_MAX) - std::log10(X_MIN)) / NUM_BINS;
    for (int i = 0; i < NUM_BINS + 1; ++i) {
        cell[i] = std::pow(10.0, std::log10(X_MIN) + i * dx);
    }
    for (int i = 0; i < NUM_BINS; ++i) {
        h[i] = cell[i + 1] - cell[i];
        x[i] = 0.5 * (cell[i + 1] + cell[i]);
    }
}

/**
 * @brief Initializes mass density distribution
 * @param g Polynomial coefficients (output)
 * @param h Cell widths
 * @param x Cell centers
 * @param cell Bin boundaries
 */
void initialize_mass_density(std::vector<std::vector<double>>& g, std::vector<double>& h, std::vector<double>& x, std::vector<double>& cell) {
    (void)cell;
    // We use the 5th order quadrature data specifically for initial projection
    const double pts5[] = {-0.9061798459, -0.5384693101, 0.0, 0.5384693101, 0.9061798459};
    const double wts5[] = {0.2369268850, 0.4786286704, 0.5688888888, 0.4786286704, 0.2369268850};
    
    for (int j = 0; j < NUM_BINS; ++j) {
        for (int i = 0; i < POLY_ORDER + 1; ++i) {
            double integral = 0.0;
            for (int k = 0; k < 5; ++k) {
                integral += wts5[k] * g_0(0.5 * h[j] * pts5[k] + x[j]) * legendre(i, pts5[k]);
            }
            g[j][i] = integral * (2.0 * i + 1.0) / 2.0;
        }
    }
    // Apply positive floor vlue to ensure non-zero mass density
    for (int j = 0; j < NUM_BINS; ++j) {
        if (std::abs(g[j][0] * h[j]) < 1e-15) {
            g[j][0] = 1e-15 * std::copysign(1.0, g[j][0]) / h[j];
            for (int i = 1; i < POLY_ORDER + 1; i++) {
                g[j][i] = 0.0;
            }
        }
    }
}

// ==========================================================================
// Functions for calculating cell averages, slope limiters, and reconstruction
// =========================================================================

/**
 * @brief Calculate cell-averaged mass density
 * @param g DG coefficients
 * @param g_bar Output: Cell averages
 */
void cal_cell_average(std::vector<std::vector<double>>& g, std::vector<double>& g_bar) {
    for (int j = 0; j < NUM_BINS; ++j) {
        g_bar[j] = g[j][0];
        if (g_bar[j] < 0) {
            std::cout << "Warning: Averaged mass density is negative at j = " << j << std::endl;
        }
    }
}

/**
 * @brief Calculates slope limiter gamma
 * @param g Polynomial coefficients
 * @param h Cell widths
 * @param x Cell centers
 * @param cell Bin boundaries
 * @param g_bar Cell averages
 * @param gamma Slope limiters (output)
 */
void cal_gamma(std::vector<std::vector<double>>& g, std::vector<double>& h, std::vector<double>& x, std::vector<double>& cell, std::vector<double>& g_bar, std::vector<double>& gamma) {
    (void)h; (void)x; (void)cell;
    cal_cell_average(g, g_bar);
    if (POLY_ORDER == 0) {
        std::fill(gamma.begin(), gamma.end(), 1.0);
        return;
    }
    for (int j = 0; j < NUM_BINS; ++j) {
        std::vector<double> extrema = {-1.0, 1.0};
        
        if (POLY_ORDER >= 2 && g[j][2] != 0) {
            double x_value = -g[j][1] / (3.0 * g[j][2]);
            if (x_value > -1.0 && x_value < 1.0) extrema.push_back(x_value);
        }
        if (POLY_ORDER >= 3 && g[j][3] != 0) {
            double a = g[j][3] * 15.0 / 2.0;
            double b = g[j][2] * 3.0;
            double c = g[j][1] - 3.0 / 2.0 * g[j][3];
            double discriminant = b * b - 4.0 * a * c;
            if (discriminant >= 0) {
                double sqrt_discriminant = std::sqrt(discriminant);
                double xi1 = (-b - sqrt_discriminant) / (2.0 * a);
                double xi2 = (-b + sqrt_discriminant) / (2.0 * a);
                if (xi1 > -1.0 && xi1 < 1.0) extrema.push_back(xi1);
                if (xi2 > -1.0 && xi2 < 1.0) extrema.push_back(xi2);
            }
        }
        if (POLY_ORDER >= 4 && g[j][4] != 0) {
            double a = g[j][4] * 35.0 / 2.0;
            double b = g[j][3] * 15.0 / 2.0;
            double c = 3.0 * g[j][2] - 15.0 / 2.0 * g[j][4];
            double d = g[j][1] - 3.0 / 2.0 * g[j][3];
            double y1 = (3.0 * a * c - b * b) / (3.0 * a * a);
            double y0 = (2.0 * b * b * b - 9.0 * a * b * c + 27.0 * a * a * d) / (27.0 * a * a * a);
            double delta = (y0 * y0 / 4.0) + (y1 * y1 * y1 / 27.0);
            if (delta > 0.0) {
                double u = std::cbrt(-y0 / 2.0 + std::sqrt(delta));
                double v = std::cbrt(-y0 / 2.0 - std::sqrt(delta));
                double root1 = u + v - b / (3.0 * a);
                if (root1 > -1.0 && root1 < 1.0) extrema.push_back(root1);
            } else if (delta == 0.0) {
                double root1 = 3.0 * y0 / y1 - b / (3.0 * a);
                double root2 = -3.0 * y0 / (2.0 * y1) - b / (3.0 * a);
                if (root1 > -1.0 && root1 < 1.0) extrema.push_back(root1);
                if (root2 > -1.0 && root2 < 1.0) extrema.push_back(root2);
            } else {
                double theta = std::acos(-y0 / 2.0 * std::sqrt(-27.0 / (y1 * y1 * y1)));
                double root1 = 2.0 * std::sqrt(-y1 / 3.0) * std::cos(theta / 3.0) - b / (3.0 * a);
                double root2 = 2.0 * std::sqrt(-y1 / 3.0) * std::cos((theta + 2.0 * M_PI) / 3.0) - b / (3.0 * a);
                double root3 = 2.0 * std::sqrt(-y1 / 3.0) * std::cos((theta + 4.0 * M_PI) / 3.0) - b / (3.0 * a);
                if (root1 > -1.0 && root1 < 1.0) extrema.push_back(root1);
                if (root2 > -1.0 && root2 < 1.0) extrema.push_back(root2);
                if (root3 > -1.0 && root3 < 1.0) extrema.push_back(root3);
            }
        }
        double g_min = std::numeric_limits<double>::max();
        for (double x_val : extrema) {
            double poly_val = 0.0;
            for (int i = 0; i <= POLY_ORDER; ++i) {
                poly_val += g[j][i] * legendre(i, x_val);
            }
            if (poly_val < g_min) g_min = poly_val;
        }
        gamma[j] = std::min(1.0, g_bar[j] / (g_bar[j] - g_min));
    }
}

/**
 * @brief Computes p: reconstruction
 * @param g Input polynomial coefficients
 * @param gamma Slope limiters
 * @param p Scaled coefficients (output)
 */
void cal_p(std::vector<std::vector<double>>& g, std::vector<double>& gamma, std::vector<std::vector<double>>& p) {
    for (int j = 0; j < NUM_BINS; ++j) {
        p[j][0] = g[j][0]; 
        for (int i = 1; i < POLY_ORDER + 1; ++i) {
            p[j][i] = g[j][i] * gamma[j];
        }
    }
}

// ============================================================================
// Flux calculations for coagulation and fragmentation
// ===========================================================================

/**
 * @brief Calculates coagulation flux at cell boundaries (non-conservative version)
 * @param p Scaled polynomial coefficients
 * @param h Cell widths
 * @param x Cell centers
 * @param cell Bin boundaries
 * @param j Boundary index
 * @return Flux value at boundary xjp1o2
 */
double cal_flux_coag_nc(std::vector<std::vector<double>>& p, std::vector<double>& h, std::vector<double>& x, std::vector<double>& cell, int j, const std::vector<std::vector<double>>& laqp) {
    if (j == -1) return 0.0;
    double flux = 0.0;
    for (int l = 0; l < j + 1; l++) {
        double flux_l = 0.0;
        for (int alpha = 0; alpha < Q_SIZE; alpha++) {
            double g_l_alpha = 0.0;
            double Gamma_l_alpha = 0.0;
            for (int i = 0; i < POLY_ORDER + 1; i++) g_l_alpha += p[l][i] * laqp[alpha][i];
            
            double x_l_alpha = x[l] + 0.5 * h[l] * Q_POINTS[alpha];
            double aj = cell[j + 1] - x_l_alpha + cell[0];
            int A = find_bin_fast(cell, aj);
            
            for (int m = A; m < NUM_BINS; m++) {
                double am = cell[m], bm = cell[m + 1];
                if (m == A) am = aj;
                double flux_m = 0.0;
                double halfplusm = 0.5 * (bm + am);
                double halfminusm = 0.5 * (bm - am);
                double invhm = 2.0 / h[m];
                
                for (int beta = 0; beta < Q_SIZE; beta++) {
                    double x_m_beta = halfplusm + halfminusm * Q_POINTS[beta];
                    double xi_m_beta = (x_m_beta - x[m]) * invhm;
                    double g_m_beta = 0.0;
                    for (int i = 0; i < POLY_ORDER + 1; i++) {
                        g_m_beta += p[m][i] * legendre(i, xi_m_beta);
                    }
                    flux_m += Q_WEIGHTS[beta] * coagkernel(x_l_alpha, x_m_beta, coagulationKernel) * g_m_beta / x_m_beta;
                }
                Gamma_l_alpha += halfminusm * flux_m;
            }
            flux_l += Q_WEIGHTS[alpha] * g_l_alpha * Gamma_l_alpha;
        }
        flux += 0.5 * h[l] * flux_l;
    }
    return flux;
}

/**
 * @brief Calculates coagulation flux integral (non-conservative version)
 * @param p Scaled polynomial coefficients
 * @param h Cell widths
 * @param x Cell centers
 * @param cell Bin boundaries
 * @param j Cell index
 * @param i Polynomial order index
 * @return Flux integral over Ij
 */
double cal_flux_integral_coag_nc(std::vector<std::vector<double>>& p, std::vector<double>& h, std::vector<double>& x, std::vector<double>& cell, int j, int i, const std::vector<std::vector<double>>& Dlaqp) {
    if (i == 0) return 0.0;
    double flux_integral = 0.0;
    for (int theta = 0;  theta < Q_SIZE; theta++) {
        double flux_j_theta = 0.0;
        double x_j_theta = x[j] + 0.5 * h[j] * Q_POINTS[theta];
        for (int l = 0; l < j + 1; l++) {
            double al = cell[l], bl = cell[l + 1];
            if (l == j) bl = x_j_theta;
            double halfplusl = 0.5 * (bl + al);
            double halfminusl = 0.5 * (bl - al);
            double invhl = 2.0 / h[l];
            double flux_l = 0.0;
            
            for (int alpha = 0; alpha < Q_SIZE; alpha++) {
                double x_l_alpha = halfplusl + halfminusl * Q_POINTS[alpha];
                double xi_l_alpha = (x_l_alpha - x[l]) * invhl;
                double g_l_alpha = 0.0, Gamma_l_alpha = 0.0;
                for (int k = 0; k < POLY_ORDER + 1; k++) g_l_alpha += p[l][k] * legendre(k, xi_l_alpha);
                
                double aj = x_j_theta - x_l_alpha + cell[0];
                int A = find_bin_fast(cell, aj);
                for (int m = A; m < NUM_BINS; m++) {
                    double am = cell[m], bm = cell[m + 1];
                    if (m == A) am = aj;
                    double flux_m = 0.0;
                    double halfplusm = 0.5 * (bm + am);
                    double halfminusm = 0.5 * (bm - am);
                    double invhm = 2.0 / h[m];
                    
                    for (int beta = 0; beta < Q_SIZE; beta++) {
                        double x_m_beta = halfplusm + halfminusm * Q_POINTS[beta];
                        double xi_m_beta = (x_m_beta - x[m]) * invhm;
                        double g_m_beta = 0.0;
                        for (int k = 0; k < POLY_ORDER + 1; k++) g_m_beta += p[m][k] * legendre(k, xi_m_beta);
                        flux_m += Q_WEIGHTS[beta] * coagkernel(x_l_alpha, x_m_beta, coagulationKernel) * g_m_beta / x_m_beta;
                    }
                    Gamma_l_alpha += halfminusm * flux_m;
                }
                flux_l += Q_WEIGHTS[alpha] * g_l_alpha * Gamma_l_alpha;
            }
            flux_j_theta += halfminusl * flux_l;
        }
        flux_integral += Q_WEIGHTS[theta] * Dlaqp[theta][i] * flux_j_theta;
    }
    return flux_integral;
}

/**
 * @brief Calculates fragmentation flux at cell boundaries (non-conservative version)
 * @param p Scaled polynomial coefficients
 * @param h Cell widths
 * @param x Cell centers
 * @param cell Bin boundaries
 * @param j Boundary index
 * @return Flux value at boundary xjp1o2
 */
double cal_flux_frag_nc(std::vector<std::vector<double> >& p, std::vector<double>& h, std::vector<double>& x, std::vector<double>& cell, int j, const std::vector<std::vector<double> >& laqp) {
    if (j == -1) return 0.0;

    double flux1 = 0.0;
    for (int l = 0; l < j + 1; l++) {
        double flux_l = 0.0;
        for (int alpha = 0; alpha < Q_SIZE; alpha++) {
            double x_l_alpha = x[l] + 0.5 * h[l] * Q_POINTS[alpha];
            double g_l_alpha = 0.0;
            for (int i = 0; i < POLY_ORDER + 1; i++) {
                g_l_alpha += p[l][i] * laqp[alpha][i];
            }
            double Gamma_l_alpha = 0.0;
            double aj = cell[j + 1] - x_l_alpha + cell[0];
            int A = find_bin_fast(cell, aj);
            for (int m = A; m < NUM_BINS; m++) {
                double flux_m = 0.0;
                double bm = cell[m + 1];
                double am = cell[m];
                if (m == A) {
                    am = aj;
                }
                double halfminusm = 0.5 * (bm - am);
                double halfplusm = 0.5 * (bm + am);
                double invhm = 2.0 / h[m];
                for (int beta = 0; beta < Q_SIZE; beta++) {
                    double x_m_beta = halfplusm + halfminusm * Q_POINTS[beta];
                    double xi_m_beta = (x_m_beta - x[m]) * invhm;
                    double g_m_beta = 0.0;
                    for (int i = 0; i < POLY_ORDER + 1; i++) {
                        g_m_beta += p[m][i] * legendre(i, xi_m_beta);
                    }
                    double integral = inte1(x_l_alpha, x_m_beta, cell[j + 1]);
                    flux_m += Q_WEIGHTS[beta] * g_m_beta / x_m_beta / (x_m_beta + x_l_alpha) * integral * fragkernel(x_l_alpha, x_m_beta, fragmentationKernel);
                }
                Gamma_l_alpha += halfminusm * flux_m;
            }
            flux_l += Q_WEIGHTS[alpha] * g_l_alpha * Gamma_l_alpha;
        }
        flux1 += h[l] * 0.5 * flux_l;
    }

    double flux2 = 0.0;
    for (int l = j + 1; l < NUM_BINS; l++) {
        double flux_l = 0.0;
        for (int alpha = 0; alpha < Q_SIZE; alpha++) {
            double x_l_alpha = x[l] + 0.5 * h[l] * Q_POINTS[alpha];
            double g_l_alpha = 0.0;
            for (int i = 0; i < POLY_ORDER + 1; i++) {
                g_l_alpha += p[l][i] * laqp[alpha][i];
            }
            double Gamma_l_alpha = 0.0;
            for (int m = 0; m < NUM_BINS; m++) {
                double flux_m = 0.0;
                for (int beta = 0; beta < Q_SIZE; beta++) {
                    double x_m_beta = x[m] + 0.5 * h[m] * Q_POINTS[beta];
                    double g_m_beta = 0.0;
                    for (int i = 0; i < POLY_ORDER + 1; i++) {
                        g_m_beta += p[m][i] * legendre(i, Q_POINTS[beta]);
                    }
                    double integral = inte2(x_l_alpha, x_m_beta, cell[j + 1]);
                    flux_m += Q_WEIGHTS[beta] * g_m_beta / x_m_beta / (x_m_beta + x_l_alpha) * integral * fragkernel(x_l_alpha, x_m_beta, fragmentationKernel);
                }
                Gamma_l_alpha += 0.5 * h[m] * flux_m;
            }
            flux_l += Q_WEIGHTS[alpha] * g_l_alpha * Gamma_l_alpha;
        }
        flux2 += h[l] * 0.5 * flux_l;
    }
    return flux1 - flux2;
}

/**
 * @brief Calculates fragmentation flux integral (non-conservative version)
 * @param p Scaled polynomial coefficients
 * @param h Cell widths
 * @param x Cell centers
 * @param cell Bin boundaries
 * @param j Cell index
 * @param i Polynomial order index
 * @return Flux integral over Ij
 */
double cal_flux_integral_frag_nc(std::vector<std::vector<double> >& p, std::vector<double>& h, std::vector<double>& x, std::vector<double>& cell, int j, int i, const std::vector<std::vector<double> >& laqp) {
    if (i == 0) return 0.0;

    double flux1 = 0.0;
    for (int l = 0; l < NUM_BINS; l++) {
        double flux_l = 0.0;
        for (int alpha = 0; alpha < Q_SIZE; alpha++) {
            double x_l_alpha = x[l] + 0.5 * h[l] * Q_POINTS[alpha];
            double g_l_alpha = 0.0;
            for (int k = 0; k < POLY_ORDER + 1; k++) {
                g_l_alpha += p[l][k] * laqp[alpha][k];
            }
            double Gamma_l_alpha = 0.0;
            for (int m = 0; m < NUM_BINS; m++) {
                double flux_m = 0.0;
                double bm = cell[m + 1];
                double am = cell[m];
                double halfplusm = 0.5 * (bm + am);
                double halfminusm = 0.5 * (bm - am);
                double invhm = 2.0 / h[m];
                for (int beta = 0; beta < Q_SIZE; beta++) {
                    double x_m_beta = halfplusm + halfminusm * Q_POINTS[beta];
                    double xi_m_beta = (x_m_beta - x[m]) * invhm;
                    double g_m_beta = 0.0;
                    for (int k = 0; k < POLY_ORDER + 1; k++) {
                        g_m_beta += p[m][k] * legendre(k, xi_m_beta);
                    }
                    double integral = 0.0;
                    double mfrag = m_frag(x_l_alpha, x_m_beta);
                    double u_v_xmin = x_l_alpha + x_m_beta - cell[0];
                    double upper = std::min(std::min(mfrag, u_v_xmin), cell[j+1]);
                    double lower = std::max(x_l_alpha, cell[j]);
                    if (lower < upper) {
                        double xi_l = (lower - x[j]) * 2.0 / h[j];
                        double xi_u = (upper - x[j]) * 2.0 / h[j];
                        double inner1 = std::pow(mfrag, alpha_p2) * (legendre(i, xi_u) - legendre(i, xi_l));
                        double inner2 = 0.0;
                        double halfminusl = 0.5 * (xi_u - xi_l);
                        double halfplusl = 0.5 * (xi_u + xi_l);
                        for (int theta = 0;  theta < Q_SIZE; theta++) {
                            double xi_j_theta = halfplusl + halfminusl * Q_POINTS[theta];
                            double x_j_theta = 0.5 * h[j] * xi_j_theta + x[j];
                            inner2 += halfminusl * Q_WEIGHTS[theta] * Dlegendre(i, xi_j_theta) * std::pow(x_j_theta, alpha_p2);
                        }
                        integral += mfrag / (std::pow(mfrag, alpha_p2) - xmin_pow) * (inner1 - inner2);
                    }
                    flux_m += Q_WEIGHTS[beta] * g_m_beta / x_m_beta / (x_m_beta + x_l_alpha) * fragkernel(x_l_alpha, x_m_beta, fragmentationKernel) * integral;
                }
                Gamma_l_alpha += halfminusm * flux_m;
            }
            flux_l += Q_WEIGHTS[alpha] * g_l_alpha * Gamma_l_alpha;
        }
        flux1 += h[l] * 0.5 * flux_l;
    }

    double flux2 = 0.0;
    for (int l = 0; l < NUM_BINS; l++) {
        double flux_l = 0.0;
        for (int alpha = 0; alpha < Q_SIZE; alpha++) {
            double x_l_alpha = x[l] + 0.5 * h[l] * Q_POINTS[alpha];
            double g_l_alpha = 0.0;
            for (int k = 0; k < POLY_ORDER + 1; k++) {
                g_l_alpha += p[l][k] * laqp[alpha][k];
            }
            double Gamma_l_alpha = 0.0;
            for (int m = 0; m < NUM_BINS; m++) {
                double flux_m = 0.0;
                for (int beta = 0; beta < Q_SIZE; beta++) {
                    double x_m_beta = x[m] + 0.5 * h[m] * Q_POINTS[beta];
                    double g_m_beta = 0.0;
                    for (int k = 0; k < POLY_ORDER + 1; k++) {
                        g_m_beta += p[m][k] * laqp[beta][k]; // Note: using laqp[beta][k] based on original main.cpp
                    }
                    double integral = 0.0;
                    double mfrag = m_frag(x_l_alpha, x_m_beta);
                    double upper = std::min(x_l_alpha, cell[j+1]);
                    double lower = cell[j];
                    if (lower < upper) {
                        double xi_l = (lower - x[j]) * 2.0 / h[j];
                        double xi_u = (upper - x[j]) * 2.0 / h[j];
                        double inner1 = xmin_pow * (legendre(i, xi_u) - legendre(i, xi_l));
                        double inner2 = 0.0;
                        double halfminusl = 0.5 * (xi_u - xi_l);
                        double halfplusl = 0.5 * (xi_u + xi_l);
                        for (int theta = 0;  theta < Q_SIZE; theta++) {
                            double xi_j_theta = halfplusl + halfminusl * Q_POINTS[theta];
                            double x_j_theta = 0.5 * h[j] * xi_j_theta + x[j];
                            inner2 += halfminusl * Q_WEIGHTS[theta] * Dlegendre(i, xi_j_theta) * std::pow(x_j_theta, alpha_p2);
                        }
                        integral += mfrag / (std::pow(mfrag, alpha_p2) - xmin_pow) * (inner2 - inner1);
                    }

                    lower = std::max(mfrag, cell[j]);
                    upper = std::min(x_l_alpha, cell[j+1]);
                    if (lower < upper) {
                        double xi_l = (lower - x[j]) * 2.0 / h[j];
                        double xi_u = (upper - x[j]) * 2.0 / h[j];
                        double inner1 = std::pow(mfrag, alpha_p2) * (legendre(i, xi_u) - legendre(i, xi_l));
                        double inner2 = 0.0;
                        double halfminusl = 0.5 * (xi_u - xi_l);
                        double halfplusl = 0.5 * (xi_u + xi_l);
                        for (int theta = 0;  theta < Q_SIZE; theta++) {
                            double xi_j_theta = halfplusl + halfminusl * Q_POINTS[theta];
                            double x_j_theta = 0.5 * h[j] * xi_j_theta + x[j];
                            inner2 += halfminusl * Q_WEIGHTS[theta] * Dlegendre(i, xi_j_theta) * std::pow(x_j_theta, alpha_p2);
                        }
                        integral += mfrag / (std::pow(mfrag, alpha_p2) - xmin_pow) * (inner1 - inner2);
                    }
                    flux_m += Q_WEIGHTS[beta] * g_m_beta / x_m_beta / (x_m_beta + x_l_alpha) * fragkernel(x_l_alpha, x_m_beta, fragmentationKernel) * integral;
                }
                Gamma_l_alpha += 0.5 * h[m] * flux_m;
            }
            flux_l += Q_WEIGHTS[alpha] * g_l_alpha * Gamma_l_alpha;
        }
        flux2 += h[l] * 0.5 * flux_l;
    }

    double flux3 = 0.0;
    for (int l = 0; l < NUM_BINS; l++) {
        double flux_l = 0.0;
        for (int alpha = 0; alpha < Q_SIZE; alpha++) {
            double x_l_alpha = x[l] + 0.5 * h[l] * Q_POINTS[alpha];
            double g_l_alpha = 0.0;
            for (int k = 0; k < POLY_ORDER + 1; k++) {
                g_l_alpha += p[l][k] * laqp[alpha][k];
            }
            double Gamma_l_alpha = 0.0;
            for (int m = 0; m < NUM_BINS; m++) {
                double flux_m = 0.0;
                double bm = cell[m + 1];
                double am = cell[m];
                double halfplusm = 0.5 * (bm + am);
                double halfminusm = 0.5 * (bm - am);
                double invhm = 2.0 / h[m];
                for (int beta = 0; beta < Q_SIZE; beta++) {
                    double x_m_beta = halfplusm + halfminusm * Q_POINTS[beta];
                    double xi_m_beta = (x_m_beta - x[m]) * invhm;
                    double g_m_beta = 0.0;
                    for (int k = 0; k < POLY_ORDER + 1; k++) {
                        g_m_beta += p[m][k] * legendre(k, xi_m_beta);
                    }
                    double integral = 0.0;
                    double mfrag = m_frag(x_l_alpha, x_m_beta);
                    double mleft = x_m_beta + x_l_alpha - mfrag;
                    double u_v_xmin = x_l_alpha + x_m_beta - cell[0];
                    double upper = std::min(std::min(mleft, u_v_xmin), cell[j+1]);
                    double lower = std::max(x_l_alpha, cell[j]);
                    if (lower < upper) {
                        double xi_l = (lower - x[j]) * 2.0 / h[j];
                        double xi_u = (upper - x[j]) * 2.0 / h[j];
                        integral += legendre(i, xi_u) - legendre(i, xi_l);
                    }
                    flux_m += Q_WEIGHTS[beta] * g_m_beta * mleft / x_m_beta / (x_m_beta + x_l_alpha) * fragkernel(x_l_alpha, x_m_beta, fragmentationKernel) * integral;
                }
                Gamma_l_alpha += halfminusm * flux_m;
            }
            flux_l += Q_WEIGHTS[alpha] * g_l_alpha * Gamma_l_alpha;
        }
        flux3 += h[l] * 0.5 * flux_l;
    }

    double flux4 = 0.0;
    for (int l = 0; l < NUM_BINS; l++) {
        double flux_l = 0.0;
        for (int alpha = 0; alpha < Q_SIZE; alpha++) {
            double x_l_alpha = x[l] + 0.5 * h[l] * Q_POINTS[alpha];
            double g_l_alpha = 0.0;
            for (int k = 0; k < POLY_ORDER + 1; k++) {
                g_l_alpha += p[l][k] * laqp[alpha][k];
            }
            double Gamma_l_alpha = 0.0;
            for (int m = 0; m < NUM_BINS; m++) {
                double flux_m = 0.0;
                for (int beta = 0; beta < Q_SIZE; beta++) {
                    double x_m_beta = x[m] + 0.5 * h[m] * Q_POINTS[beta];
                    double g_m_beta = 0.0;
                    for (int k = 0; k < POLY_ORDER + 1; k++) {
                        g_m_beta += p[m][k] * laqp[beta][k];
                    }
                    double integral = 0.0;
                    double mfrag = m_frag(x_l_alpha, x_m_beta);
                    double mleft = x_m_beta + x_l_alpha - mfrag;
                    double lower = std::max(mleft, cell[j]);
                    double upper = std::min(x_l_alpha, cell[j+1]);
                    if (lower < upper) {
                        double xi_l = (lower - x[j]) * 2.0 / h[j];
                        double xi_u = (upper - x[j]) * 2.0 / h[j];
                        integral += legendre(i, xi_u) - legendre(i, xi_l);
                    }
                    flux_m += Q_WEIGHTS[beta] * g_m_beta * mleft / x_m_beta / (x_m_beta + x_l_alpha) * fragkernel(x_l_alpha, x_m_beta, fragmentationKernel) * integral;
                }
                Gamma_l_alpha += 0.5 * h[m] * flux_m;
            }
            flux_l += Q_WEIGHTS[alpha] * g_l_alpha * Gamma_l_alpha;
        }
        flux4 += h[l] * 0.5 * flux_l;
    }
    double flux_integral = flux1 - flux2 + flux3 - flux4;
    return flux_integral;
}

/**
 * @brief Calculates coagulation flux at cell boundaries (conservative version)
 * @param p Scaled polynomial coefficients
 * @param h Cell widths
 * @param x Cell centers
 * @param cell Bin boundaries
 * @param j Boundary index
 * @return Flux value at boundary xjp1o2
 */
double cal_flux_coag_c(std::vector<std::vector<double> >& p, std::vector<double>& h, std::vector<double>& x, std::vector<double>& cell, int j, const std::vector<std::vector<double> >& laqp) {
    if (j == -1) return 0.0;
    if (j == NUM_BINS - 1) return 0.0;

    double flux = 0.0;
    for (int l = 0; l < j + 1; l++) {
        double flux_l = 0.0;
        for (int alpha = 0; alpha < Q_SIZE; alpha++) {
            double x_l_alpha = x[l] + 0.5 * h[l] * Q_POINTS[alpha];
            double g_l_alpha = 0.0;
            double Gamma_l_alpha = 0.0;
            for (int i = 0; i < POLY_ORDER + 1; i++) {
                g_l_alpha += p[l][i] * laqp[alpha][i];
            }
            double bj = cell[NUM_BINS] - x_l_alpha + cell[0];
            double aj = cell[j + 1] - x_l_alpha + cell[0];
            int B = find_bin_fast(cell, bj);
            int A = find_bin_fast(cell, aj);
            for (int m = A; m < B + 1; m++) {
                double am = cell[m];
                double bm = cell[m + 1];
                if (m == A) {
                    am = aj;
                }
                if (m == B) {
                    bm = bj;
                }
                double flux_m = 0.0;
                double halfplusm = 0.5 * (bm + am);
                double halfminusm = 0.5 * (bm - am);
                double invhm = 2.0 / h[m];
                for (int beta = 0; beta < Q_SIZE; beta++) {
                    double x_m_beta = halfplusm + halfminusm * Q_POINTS[beta];
                    double xi_m_beta = (x_m_beta - x[m]) * invhm;
                    double g_m_beta = 0.0;
                    for (int i = 0; i < POLY_ORDER + 1; i++) {
                        g_m_beta += p[m][i] * legendre(i, xi_m_beta);
                    }
                    flux_m += Q_WEIGHTS[beta] * coagkernel(x_l_alpha, x_m_beta, coagulationKernel) * g_m_beta / x_m_beta;
                }
                Gamma_l_alpha += halfminusm * flux_m;
            }
            flux_l += Q_WEIGHTS[alpha] * g_l_alpha * Gamma_l_alpha;
        }
        flux += 0.5 * h[l] * flux_l;
    }
    return flux;
}

/**
 * @brief Calculates coagulation flux integral (conservative version)
 * @param p Scaled polynomial coefficients
 * @param h Cell widths
 * @param x Cell centers
 * @param cell Bin boundaries
 * @param j Cell index
 * @param i Polynomial order index
 * @return Flux integral over Ij
 */
double cal_flux_integral_coag_c(std::vector<std::vector<double> >& p, std::vector<double>& h, std::vector<double>& x, std::vector<double>& cell, int j, int i, const std::vector<std::vector<double> >& Dlaqp) {
    if (i == 0) return 0.0;

    double flux_integral = 0.0;
    for (int theta = 0;  theta < Q_SIZE; theta++) {
        double flux_j_theta = 0.0;
        double x_j_theta = x[j] + 0.5 * h[j] * Q_POINTS[theta];
        for (int l = 0; l < j + 1; l++) {
            double al = cell[l];
            double bl = cell[l + 1];
            if (l == j) {
                bl = x_j_theta;
            }
            double flux_l = 0.0;
            double halfplusl = 0.5 * (bl + al);
            double halfminusl = 0.5 * (bl - al);
            double invhl = 2.0 / h[l];
            for (int alpha = 0; alpha < Q_SIZE; alpha++) {
                double x_l_alpha = halfplusl + halfminusl * Q_POINTS[alpha];
                double xi_l_alpha = (x_l_alpha - x[l]) * invhl;
                double g_l_alpha = 0.0;
                double Gamma_l_alpha = 0.0;
                for (int k = 0; k < POLY_ORDER + 1; k++) {
                    g_l_alpha += p[l][k] * legendre(k, xi_l_alpha);
                }
                double bj = cell[NUM_BINS] - x_l_alpha + cell[0];
                double aj = x_j_theta - x_l_alpha + cell[0];
                int B = find_bin_fast(cell, bj);
                int A = find_bin_fast(cell, aj);
                for (int m = A; m < B + 1; m++) {
                    double am = cell[m];
                    double bm = cell[m + 1];
                    if (m == A) {
                        am = aj;
                    }
                    if (m == B) {
                        bm = bj;
                    }
                    double flux_m = 0.0;
                    double halfplusm = 0.5 * (bm + am);
                    double halfminusm = 0.5 * (bm - am);
                    double invhm = 2.0 / h[m];
                    for (int beta = 0; beta < Q_SIZE; beta++) {
                        double x_m_beta = halfplusm + halfminusm * Q_POINTS[beta];
                        double xi_m_beta = (x_m_beta - x[m]) * invhm;
                        double g_m_beta = 0.0;
                        for (int k = 0; k < POLY_ORDER + 1; k++) {
                            g_m_beta += p[m][k] * legendre(k, xi_m_beta);
                        }
                        flux_m += Q_WEIGHTS[beta] * coagkernel(x_l_alpha, x_m_beta, coagulationKernel) * g_m_beta / x_m_beta;
                    }
                    Gamma_l_alpha += halfminusm * flux_m;
                }
                flux_l += Q_WEIGHTS[alpha] * g_l_alpha * Gamma_l_alpha;
            }
            flux_j_theta += halfminusl * flux_l;
        }
        flux_integral += Q_WEIGHTS[theta] * Dlaqp[theta][i] * flux_j_theta;
    }
    return flux_integral;
}

/**
 * @brief Calculates fragmentation flux at cell boundaries (conservative version)
 * @param p Scaled polynomial coefficients
 * @param h Cell widths
 * @param x Cell centers
 * @param cell Bin boundaries
 * @param j Boundary index
 * @return Flux value at boundary xjp1o2
 */
double cal_flux_frag_c(std::vector<std::vector<double> >& p, std::vector<double>& h, std::vector<double>& x, std::vector<double>& cell, int j, const std::vector<std::vector<double> >& laqp) {
    if (j == -1) return 0.0;
    if (j == NUM_BINS - 1) return 0.0;

    double flux1 = 0.0;
    for (int l = 0; l < j + 1; l++) {
        double flux_l = 0.0;
        for (int alpha = 0; alpha < Q_SIZE; alpha++) {
            double x_l_alpha = x[l] + 0.5 * h[l] * Q_POINTS[alpha];
            double g_l_alpha = 0.0;
            for (int i = 0; i < POLY_ORDER + 1; i++) {
                g_l_alpha += p[l][i] * laqp[alpha][i];
            }
            double Gamma_l_alpha = 0.0;
            double bj = cell[NUM_BINS] - x_l_alpha + cell[0];
            double aj = cell[j + 1] - x_l_alpha + cell[0];
            int B = find_bin_fast(cell, bj);
            int A = find_bin_fast(cell, aj);
            for (int m = A; m < B + 1; m++) {
                double bm = cell[m + 1];
                double am = cell[m];
                if (m == A) {
                    am = aj;
                }
                if (m == B) {
                    bm = bj;
                }
                double flux_m = 0.0;
                double halfplusm = 0.5 * (bm + am);
                double halfminusm = 0.5 * (bm - am);
                double invhm = 2.0 / h[m];
                for (int beta = 0; beta < Q_SIZE; beta++) {
                    double x_m_beta = halfplusm + halfminusm * Q_POINTS[beta];
                    double xi_m_beta = (x_m_beta - x[m]) * invhm;
                    double g_m_beta = 0.0;
                    for (int i = 0; i < POLY_ORDER + 1; i++) {
                        g_m_beta += p[m][i] * legendre(i, xi_m_beta);
                    }
                    double integral = inte1(x_l_alpha, x_m_beta, cell[j + 1]);
                    flux_m += Q_WEIGHTS[beta] * g_m_beta / x_m_beta / (x_m_beta + x_l_alpha) * integral * fragkernel(x_l_alpha, x_m_beta, fragmentationKernel);
                }
                Gamma_l_alpha += halfminusm * flux_m;
            }
            flux_l += Q_WEIGHTS[alpha] * g_l_alpha * Gamma_l_alpha;
        }
        flux1 += h[l] * 0.5 * flux_l;
    }

    double flux2 = 0.0;
    for (int l = j + 1; l < NUM_BINS; l++) {
        double flux_l = 0.0;
        for (int alpha = 0; alpha < Q_SIZE; alpha++) {
            double x_l_alpha = x[l] + 0.5 * h[l] * Q_POINTS[alpha];
            double g_l_alpha = 0.0;
            for (int i = 0; i < POLY_ORDER + 1; i++) {
                g_l_alpha += p[l][i] * laqp[alpha][i];
            }
            double Gamma_l_alpha = 0.0;
            double bj = cell[NUM_BINS] - x_l_alpha + cell[0];
            int B = find_bin_fast(cell, bj);
            for (int m = 0; m < B + 1; m++) {
                double bm = cell[m + 1];
                double am = cell[m];
                if (m == B) {
                    bm = bj;
                }
                double flux_m = 0.0;
                double halfplusm = 0.5 * (bm + am);
                double halfminusm = 0.5 * (bm - am);
                double invhm = 2.0 / h[m];
                for (int beta = 0; beta < Q_SIZE; beta++) {
                    double x_m_beta = halfplusm + halfminusm * Q_POINTS[beta];
                    double xi_m_beta = (x_m_beta - x[m]) * invhm;
                    double g_m_beta = 0.0;
                    for (int i = 0; i < POLY_ORDER + 1; i++) {
                        g_m_beta += p[m][i] * legendre(i, xi_m_beta);
                    }
                    double integral = inte2(x_l_alpha, x_m_beta, cell[j + 1]);
                    flux_m += Q_WEIGHTS[beta] * g_m_beta / x_m_beta / (x_m_beta + x_l_alpha) * integral * fragkernel(x_l_alpha, x_m_beta, fragmentationKernel);
                }
                Gamma_l_alpha += halfminusm * flux_m;
            }
            flux_l += Q_WEIGHTS[alpha] * g_l_alpha * Gamma_l_alpha;
        }
        flux2 += h[l] * 0.5 * flux_l;
    }
    return flux1 - flux2;
}

/**
 * @brief Calculates fragmentation flux integral (conservative version)
 * @param p Scaled polynomial coefficients
 * @param h Cell widths
 * @param x Cell centers
 * @param cell Bin boundaries
 * @param j Cell index
 * @param i Polynomial order index
 * @return Flux integral over Ij
 */
double cal_flux_integral_frag_c(std::vector<std::vector<double> >& p, std::vector<double>& h, std::vector<double>& x, std::vector<double>& cell, int j, int i, const std::vector<std::vector<double> >& laqp) {
    if (i == 0) return 0.0;

    double flux1 = 0.0;
    for (int l = 0; l < NUM_BINS; l++) {
        double flux_l = 0.0;
        for (int alpha = 0; alpha < Q_SIZE; alpha++) {
            double x_l_alpha = x[l] + 0.5 * h[l] * Q_POINTS[alpha];
            double g_l_alpha = 0.0;
            for (int k = 0; k < POLY_ORDER + 1; k++) {
                g_l_alpha += p[l][k] * laqp[alpha][i]; // Note: using laqp[alpha][i] from original
            }
            double Gamma_l_alpha = 0.0;
            double bj = cell[NUM_BINS] - x_l_alpha + cell[0];
            int B = find_bin_fast(cell, bj);
            for (int m = 0; m < B + 1; m++) {
                double bm = cell[m + 1];
                double am = cell[m];
                if (m == B) {
                    bm = bj;
                }
                double flux_m = 0.0;
                double halfplusm = 0.5 * (bm + am);
                double halfminusm = 0.5 * (bm - am);
                double invhm = 2.0 / h[m];
                for (int beta = 0; beta < Q_SIZE; beta++) {
                    double x_m_beta = halfplusm + halfminusm * Q_POINTS[beta];
                    double xi_m_beta = (x_m_beta - x[m]) * invhm;
                    double g_m_beta = 0.0;
                    for (int k = 0; k < POLY_ORDER + 1; k++) {
                        g_m_beta += p[m][k] * legendre(k, xi_m_beta);
                    }
                    double integral = 0.0;
                    double mfrag = m_frag(x_m_beta, x_l_alpha);
                    double u_v_xmin = x_l_alpha + x_m_beta - cell[0];
                    double upper = std::min(std::min(mfrag, u_v_xmin), cell[j+1]);
                    double lower = std::max(x_l_alpha, cell[j]);
                    if (lower < upper) {
                        double xi_l = (lower - x[j]) * 2.0 / h[j];
                        double xi_u = (upper - x[j]) * 2.0 / h[j];
                        double inner1 = std::pow(mfrag, alpha_p2) * (legendre(i, xi_u) - legendre(i, xi_l));
                        double inner2 = 0.0;
                        double halfminusl = 0.5 * (xi_u - xi_l);
                        double halfplusl = 0.5 * (xi_u + xi_l);
                        for (int theta = 0;  theta < Q_SIZE; theta++) {
                            double xi_j_theta = halfplusl + halfminusl * Q_POINTS[theta];
                            double x_j_theta = 0.5 * h[j] * xi_j_theta + x[j];
                            inner2 += halfminusl * Q_WEIGHTS[theta] * Dlegendre(i, xi_j_theta) * std::pow(x_j_theta, alpha_p2);
                        }
                        integral += mfrag / (std::pow(mfrag, alpha_p2) - xmin_pow) * (inner1 - inner2);
                    }
                    flux_m += Q_WEIGHTS[beta] * g_m_beta / x_m_beta / (x_m_beta + x_l_alpha) * fragkernel(x_l_alpha, x_m_beta, fragmentationKernel) * integral;
                }
                Gamma_l_alpha += halfminusm * flux_m;
            }
            flux_l += Q_WEIGHTS[alpha] * g_l_alpha * Gamma_l_alpha;
        }
        flux1 += h[l] * 0.5 * flux_l;
    }

    double flux2 = 0.0;
    for (int l = 0; l < NUM_BINS; l++) {
        double flux_l = 0.0;
        for (int alpha = 0; alpha < Q_SIZE; alpha++) {
            double x_l_alpha = x[l] + 0.5 * h[l] * Q_POINTS[alpha];
            double g_l_alpha = 0.0;
            for (int k = 0; k < POLY_ORDER + 1; k++) {
                g_l_alpha += p[l][k] * laqp[alpha][k];
            }
            double Gamma_l_alpha = 0.0;
            double bj = cell[NUM_BINS] - x_l_alpha + cell[0];
            int B = find_bin_fast(cell, bj);
            for (int m = 0; m < B + 1; m++) {
                double flux_m = 0.0;
                double bm = cell[m + 1];
                double am = cell[m];
                if (m == B) {
                    bm = bj;
                }
                double halfplusm = 0.5 * (bm + am);
                double halfminusm = 0.5 * (bm - am);
                double invhm = 2.0 / h[m];
                for (int beta = 0; beta < Q_SIZE; beta++) {
                    double x_m_beta = halfplusm + halfminusm * Q_POINTS[beta];
                    double xi_m_beta = (x_m_beta - x[m]) * invhm;
                    double g_m_beta = 0.0;
                    for (int k = 0; k < POLY_ORDER + 1; k++) {
                        g_m_beta += p[m][k] * legendre(k, xi_m_beta);
                    }
                    double integral = 0.0;
                    double mfrag = m_frag(x_m_beta, x_l_alpha);
                    double upper = std::min(x_l_alpha, cell[j+1]);
                    double lower = cell[j];
                    if (lower < upper) {
                        double xi_l = (lower - x[j]) * 2.0 / h[j];
                        double xi_u = (upper - x[j]) * 2.0 / h[j];
                        double inner1 = xmin_pow * (legendre(i, xi_u) - legendre(i, xi_l));
                        double inner2 = 0.0;
                        double halfminusl = 0.5 * (xi_u - xi_l);
                        double halfplusl = 0.5 * (xi_u + xi_l);
                        for (int theta = 0;  theta < Q_SIZE; theta++) {
                            double xi_j_theta = halfplusl + halfminusl * Q_POINTS[theta];
                            double x_j_theta = 0.5 * h[j] * xi_j_theta + x[j];
                            inner2 += halfminusl * Q_WEIGHTS[theta] * Dlegendre(i, xi_j_theta) * std::pow(x_j_theta, alpha_p2);
                        }
                        integral += mfrag / (std::pow(mfrag, alpha_p2) - xmin_pow) * (inner2 - inner1);
                    }

                    lower = std::max(mfrag, cell[j]);
                    upper = std::min(x_l_alpha, cell[j+1]);
                    if (lower < upper) {
                        double xi_l = (lower - x[j]) * 2.0 / h[j];
                        double xi_u = (upper - x[j]) * 2.0 / h[j];
                        double inner1 = std::pow(mfrag, alpha_p2) * (legendre(i, xi_u) - legendre(i, xi_l));
                        double inner2 = 0.0;
                        double halfminusl = 0.5 * (xi_u - xi_l);
                        double halfplusl = 0.5 * (xi_u + xi_l);
                        for (int theta = 0;  theta < Q_SIZE; theta++) {
                            double xi_j_theta = halfplusl + halfminusl * Q_POINTS[theta];
                            double x_j_theta = 0.5 * h[j] * xi_j_theta + x[j];
                            inner2 += halfminusl * Q_WEIGHTS[theta] * Dlegendre(i, xi_j_theta) * std::pow(x_j_theta, alpha_p2);
                        }
                        integral += mfrag / (std::pow(mfrag, alpha_p2) - xmin_pow) * (inner1 - inner2);
                    }
                    flux_m += Q_WEIGHTS[beta] * g_m_beta / x_m_beta / (x_m_beta + x_l_alpha) * fragkernel(x_l_alpha, x_m_beta, fragmentationKernel) * integral;
                }
                Gamma_l_alpha += halfminusm * flux_m;
            }
            flux_l += Q_WEIGHTS[alpha] * g_l_alpha * Gamma_l_alpha;
        }
        flux2 += h[l] * 0.5 * flux_l;
    }

    double flux3 = 0.0;
    for (int l = 0; l < NUM_BINS; l++) {
        double flux_l = 0.0;
        for (int alpha = 0; alpha < Q_SIZE; alpha++) {
            double x_l_alpha = x[l] + 0.5 * h[l] * Q_POINTS[alpha];
            double g_l_alpha = 0.0;
            for (int k = 0; k < POLY_ORDER + 1; k++) {
                g_l_alpha += p[l][k] * laqp[alpha][k];
            }
            double Gamma_l_alpha = 0.0;
            double bj = cell[NUM_BINS] - x_l_alpha + cell[0];
            int B = find_bin_fast(cell, bj);
            for (int m = 0; m < B + 1; m++) {
                double flux_m = 0.0;
                double bm = cell[m + 1];
                double am = cell[m];
                if (m == B) {
                    bm = bj;
                }
                double halfplusm = 0.5 * (bm + am);
                double halfminusm = 0.5 * (bm - am);
                double invhm = 2.0 / h[m];
                for (int beta = 0; beta < Q_SIZE; beta++) {
                    double x_m_beta = halfplusm + halfminusm * Q_POINTS[beta];
                    double xi_m_beta = (x_m_beta - x[m]) * invhm;
                    double g_m_beta = 0.0;
                    for (int k = 0; k < POLY_ORDER + 1; k++) {
                        g_m_beta += p[m][k] * legendre(k, xi_m_beta);
                    }
                    double integral = 0.0;
                    double mfrag = m_frag(x_m_beta, x_l_alpha);
                    double mleft = x_m_beta + x_l_alpha - mfrag;
                    double u_v_xmin = x_l_alpha + x_m_beta - cell[0];
                    double upper = std::min(std::min(mleft, u_v_xmin), cell[j+1]);
                    double lower = std::max(x_l_alpha, cell[j]);
                    if (lower < upper) {
                        double xi_l = (lower - x[j]) * 2.0 / h[j];
                        double xi_u = (upper - x[j]) * 2.0 / h[j];
                        integral += legendre(i, xi_u) - legendre(i, xi_l);
                    }
                    flux_m += Q_WEIGHTS[beta] * g_m_beta * mleft / x_m_beta / (x_m_beta + x_l_alpha) * fragkernel(x_l_alpha, x_m_beta, fragmentationKernel) * integral;
                }
                Gamma_l_alpha += halfminusm * flux_m;
            }
            flux_l += Q_WEIGHTS[alpha] * g_l_alpha * Gamma_l_alpha;
        }
        flux3 += h[l] * 0.5 * flux_l;
    }

    double flux4 = 0.0;
    for (int l = 0; l < NUM_BINS; l++) {
        double flux_l = 0.0;
        for (int alpha = 0; alpha < Q_SIZE; alpha++) {
            double x_l_alpha = x[l] + 0.5 * h[l] * Q_POINTS[alpha];
            double g_l_alpha = 0.0;
            for (int k = 0; k < POLY_ORDER + 1; k++) {
                g_l_alpha += p[l][k] * laqp[alpha][k];
            }
            double Gamma_l_alpha = 0.0;
            double bj = cell[NUM_BINS] - x_l_alpha + cell[0];
            int B = find_bin_fast(cell, bj);
            for (int m = 0; m < B + 1; m++) {
                double flux_m = 0.0;
                double bm = cell[m + 1];
                double am = cell[m];
                if (m == B) {
                    bm = bj;
                }
                double halfplusm = 0.5 * (bm + am);
                double halfminusm = 0.5 * (bm - am);
                double invhm = 2.0 / h[m];
                for (int beta = 0; beta < Q_SIZE; beta++) {
                    double x_m_beta = halfplusm + halfminusm * Q_POINTS[beta];
                    double xi_m_beta = (x_m_beta - x[m]) * invhm;
                    double g_m_beta = 0.0;
                    for (int k = 0; k < POLY_ORDER + 1; k++) {
                        g_m_beta += p[m][k] * legendre(k, xi_m_beta);
                    }
                    double integral = 0.0;
                    double mfrag = m_frag(x_m_beta, x_l_alpha);
                    double mleft = x_m_beta + x_l_alpha - mfrag;
                    double lower = std::max(mleft, cell[j]);
                    double upper = std::min(x_l_alpha, cell[j+1]);
                    if (lower < upper) {
                        double xi_l = (lower - x[j]) * 2.0 / h[j];
                        double xi_u = (upper - x[j]) * 2.0 / h[j];
                        integral += legendre(i, xi_u) - legendre(i, xi_l);
                    }
                    flux_m += Q_WEIGHTS[beta] * g_m_beta * mleft / x_m_beta / (x_m_beta + x_l_alpha) * fragkernel(x_l_alpha, x_m_beta, fragmentationKernel) * integral;
                }
                Gamma_l_alpha += halfminusm * flux_m;
            }
            flux_l += Q_WEIGHTS[alpha] * g_l_alpha * Gamma_l_alpha;
        }
        flux4 += h[l] * 0.5 * flux_l;
    }
    double flux_integral = flux1 - flux2 + flux3 - flux4;
    return flux_integral;
}

// -----------------------------------------------------------------------------
// Time-stepping and evolution Functions
// -----------------------------------------------------------------------------

/**
 * @brief Calculates adaptive time step based on flux conditions
 * @param p Scaled polynomial coefficients
 * @param h Cell widths
 * @param x Cell centers
 * @param cell Bin boundaries
 * @param flux Flux values (output)
 * @return Calculated time step
 */
double cal_time_step(std::vector<std::vector<double>>& p, std::vector<double>& h, std::vector<double>& x, std::vector<double>& cell, std::vector<double>& flux, const std::vector<std::vector<double>>& laqp) {
    std::vector<double> delta_t(NUM_BINS);
    for (int j = 0; j < NUM_BINS; j++) {
        if (flux_mode == "coag_frag") {
            flux[j] = cal_flux_coag_nc(p, h, x, cell, j-1, laqp) + cal_flux_frag_nc(p, h, x, cell, j-1, laqp);
        } else if (flux_mode == "pure_frag") {
            flux[j] = cal_flux_frag_nc(p, h, x, cell, j-1, laqp);
        } else if (flux_mode == "pure_coag") {
            flux[j] = cal_flux_coag_nc(p, h, x, cell, j-1, laqp);
        }
        
        if (j == NUM_BINS-1) {
            if (flux_mode == "coag_frag") flux[j+1] = cal_flux_coag_nc(p, h, x, cell, j, laqp) + cal_flux_frag_nc(p, h, x, cell, j, laqp);
            else if (flux_mode == "pure_frag") flux[j+1] = cal_flux_frag_nc(p, h, x, cell, j, laqp);
            else if (flux_mode == "pure_coag") flux[j+1] = cal_flux_coag_nc(p, h, x, cell, j, laqp);
        }
    }
    for (int j = 0; j < NUM_BINS; j++) {
        double flux_diff = flux[j] - flux[j+1];
        delta_t[j] = h[j] * p[j][0] / std::abs(flux_diff);
    }
    return *std::min_element(delta_t.begin(), delta_t.end()) * CFL;
}

/**
 * @brief Main evolution function
 * @param g Updated polynomial coefficients (output)
 * @param p Current polynomial coefficients
 * @param h Cell widths
 * @param x Cell centers
 * @param cell Bin boundaries
 * @param flux Flux values
 * @param dt Time step (adjusted if needed)
 * @return True if evolution successful, false if dt needs reduction
 */
bool evolve(std::vector<std::vector<double>>& g, std::vector<std::vector<double>>& p, std::vector<double>& h, std::vector<double>& x, std::vector<double>& cell, std::vector<double>& flux, double& dt, const std::vector<std::vector<double>>& laqp, const std::vector<std::vector<double>>& Dlaqp) {
    // Benckmark tests use RK3
    std::vector<std::vector<double>> g1(NUM_BINS, std::vector<double>(POLY_ORDER + 1));
    std::vector<std::vector<double>> g2(NUM_BINS, std::vector<double>(POLY_ORDER + 1));

    for (int j = 0; j < NUM_BINS; j++) {
        for (int i = 0; i < POLY_ORDER + 1; i++) {
            double update = 0.0;
            if (flux_mode == "coag_frag") update = cal_flux_integral_coag_nc(p, h, x, cell, j, i, Dlaqp) + cal_flux_integral_frag_nc(p, h, x, cell, j, i, laqp);
            else if (flux_mode == "pure_frag") update = cal_flux_integral_frag_nc(p, h, x, cell, j, i, laqp);
            else if (flux_mode == "pure_coag") update = cal_flux_integral_coag_nc(p, h, x, cell, j, i, Dlaqp);
            
            update += -flux[j+1] + flux[j] * std::pow(-1, i);
            g1[j][i] = p[j][i] + dt * (2.0 * i + 1.0) / h[j] * update;
        }
    }
    
    cal_time_step(g1, h, x, cell, flux, laqp);
    for (int j = 0; j < NUM_BINS; j++) {
        for (int i = 0; i < POLY_ORDER + 1; i++) {
            double update = 0.0;
            if (flux_mode == "coag_frag") update = cal_flux_integral_coag_nc(g1, h, x, cell, j, i, Dlaqp) + cal_flux_integral_frag_nc(g1, h, x, cell, j, i, laqp);
            else if (flux_mode == "pure_frag") update = cal_flux_integral_frag_nc(g1, h, x, cell, j, i, laqp);
            else if (flux_mode == "pure_coag") update = cal_flux_integral_coag_nc(g1, h, x, cell, j, i, Dlaqp);
            
            update += -flux[j+1] + flux[j] * std::pow(-1, i);
            g2[j][i] = 3.0 / 4.0 * p[j][i] + 1.0 / 4.0 * (g1[j][i] + dt * (2.0 * i + 1.0) / h[j] * update);
        }
    }

    cal_time_step(g2, h, x, cell, flux, laqp);
    for (int j = 0; j < NUM_BINS; j++) {
        for (int i = 0; i < POLY_ORDER + 1; i++) {
            double update = 0.0;
            if (flux_mode == "coag_frag") update = cal_flux_integral_coag_nc(g2, h, x, cell, j, i, Dlaqp) + cal_flux_integral_frag_nc(g2, h, x, cell, j, i, laqp);
            else if (flux_mode == "pure_frag") update = cal_flux_integral_frag_nc(g2, h, x, cell, j, i, laqp);
            else if (flux_mode == "pure_coag") update = cal_flux_integral_coag_nc(g2, h, x, cell, j, i, Dlaqp);
            
            update += -flux[j+1] + flux[j] * std::pow(-1, i); 
            g[j][i] = 1.0 / 3.0 * p[j][i] + 2.0 / 3.0 * (g2[j][i] + dt * (2.0 * i + 1.0) / h[j] * update);
        }
    }
    
    for (int j = 0; j < NUM_BINS; j++) {
        // Apply floor value
        if (g[j][0] * h[j] < 1e-30 && g[j][0] > 0) {
            g[j][0] = 1e-30 / h[j];
            for (int i = 1; i < POLY_ORDER + 1; i++) g[j][i] = 0.0;
        }
        // Check for negative mass densities
        if (g[j][0] < 0) {
            std::cerr << "Negative density at j=" << j << ": " << g[j][0] << std::endl;
            dt /= 2.0;
            return false;
        }
    }
    return true;
}
