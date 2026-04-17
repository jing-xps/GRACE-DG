/// This script is for the numerical benchmark tests in the code paper.
/// It implements the DG method for solving the Smoluchowski coagulation-fragmentation equation
/// with various coagulation and fragmentation kernels to validate the accuracy and efficiency of the method.
/// Author: Jing Yang

#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <fstream>
#include <math.h>
#include <chrono>

/**
 * @struct TimeProfile
 * @brief Structure for tracking computation time of different components
 */
struct TimeProfile {
    double total_cal_time_step = 0;  ///< Time for calculating time steps
    double total_evolve = 0;         ///< Time for evolution calculations
    double total_cal_gamma = 0;      ///< Time for gamma calculations
    double total_cal_p = 0;          ///< Time for reconstruction
    double total_record = 0;         ///< Time for data recording
    double total_loop = 0;           ///< Total loop time
};

// Params 
int NUM_BINS;      ///< Number of bins (cells)
int POLY_ORDER;    ///< Polynomial order (max 4)
int QUAD;          ///< Quadrature points count
double MAX_TIME;   ///< Maximum simulation time
double X_MIN;      ///< Minimum mass value
double X_MAX;      ///< Maximum mass value
double CFL;        ///< CFL number
double VELO_F;     ///< Fragmentation velocity
double VELO_B;     ///< Bouncing velocity
double ALPHA;      ///< Fragment distribution power law index

// Coagulation kernel types
enum CoagulationKernel {
    COAG_CONSTANT,       ///< Constant kernel: K = 1
    COAG_ADDITIVE,       ///< Additive kernel: K = u + v
    COAG_MULTIPLICATIVE, ///< Multiplicative kernel: K = u * v
    COAG_BALLISTIC,      ///< Ballistic kernel: K = π*(u^{1/3} + v^{1/3})^2
    COAG_DETER,          ///< Deterministic kernel with velocity thresholds
    COAG_PROB,           ///< Probabilistic kernel with velocity thresholds
};
CoagulationKernel coagulationKernel;  ///< Selected coagulation kernel

// Fragmentation kernel types (same options as coagulation)
enum FragmentationKernel {
    FRAG_CONSTANT,       ///< Constant fragmentation kernel
    FRAG_ADDITIVE,       ///< Additive fragmentation kernel
    FRAG_MULTIPLICATIVE, ///< Multiplicative fragmentation kernel
    FRAG_BALLISTIC,      ///< Ballistic fragmentation kernel
    FRAG_DETER,          ///< Deterministic fragmentation with thresholds
    FRAG_PROB,           ///< Probabilistic fragmentation with thresholds
};
FragmentationKernel fragmentationKernel;  ///< Selected fragmentation kernel

// Flux calculation modes
std::string flux_mode = "coag_frag";  ///< Default: both coagulation and fragmentation

/**
 * @namespace GaussQuadrature
 * @brief Contains Gauss quadrature points and weights for different orders
 */
namespace GaussQuadrature {
    // QUAD=1
    const double pts1[] = {0.0};
    const double wts1[] = {2.0};
    // QUAD=2
    const double pts2[] = {-1.0 / std::sqrt(3.0), 1.0 / std::sqrt(3.0)};
    const double wts2[] = {1.0, 1.0};
    // QUAD=3
    const double pts3[] = {
        -std::sqrt(3.0 / 5.0),
        0.0,
        std::sqrt(3.0 / 5.0)
    };
    const double wts3[] = {5.0 / 9.0, 8.0 / 9.0, 5.0 / 9.0};
    // QUAD=4
    const double sqrt30 = std::sqrt(30.0);
    const double pts4[] = {
        -std::sqrt((3.0 + 2.0 * sqrt30 / 5.0) / 7.0),
        -std::sqrt((3.0 - 2.0 * sqrt30 / 5.0) / 7.0),
        std::sqrt((3.0 - 2.0 * sqrt30 / 5.0) / 7.0),
        std::sqrt((3.0 + 2.0 * sqrt30 / 5.0) / 7.0)
    };
    const double wts4[] = {
        (18.0 - sqrt30) / 36.0,
        (18.0 + sqrt30) / 36.0,
        (18.0 + sqrt30) / 36.0,
        (18.0 - sqrt30) / 36.0
    };
    // QUAD=5
    const double sqrt70 = std::sqrt(70.0);
    const double pts5[] = {
        -std::sqrt(5.0 + 2.0 * sqrt70 / 7.0) / 3.0,
        -std::sqrt(5.0 - 2.0 * sqrt70 / 7.0) / 3.0,
        0.0,
        std::sqrt(5.0 - 2.0 * sqrt70 / 7.0) / 3.0,
        std::sqrt(5.0 + 2.0 * sqrt70 / 7.0) / 3.0
    };
    const double wts5[] = {
        (322.0 - 13.0 * sqrt70) / 900.0,
        (322.0 + 13.0 * sqrt70) / 900.0,
        128.0 / 225.0,
        (322.0 + 13.0 * sqrt70) / 900.0,
        (322.0 - 13.0 * sqrt70) / 900.0
    };
}

/**
 * @struct QuadratureData
 * @brief Stores quadrature points, weights and size
 */
struct QuadratureData {
    const double* points;   ///< Quadrature points array
    const double* weights;  ///< Quadrature weights array
    int quad_size;          ///< Number of quadrature points
};

/// Quadrature data for orders 1-5
const QuadratureData quadratureData[] = {
    {nullptr, nullptr, 0}, // QUAD=0 (invalid)
    {GaussQuadrature::pts1, GaussQuadrature::wts1, 1}, // QUAD=1
    {GaussQuadrature::pts2, GaussQuadrature::wts2, 2}, // QUAD=2
    {GaussQuadrature::pts3, GaussQuadrature::wts3, 3}, // QUAD=3
    {GaussQuadrature::pts4, GaussQuadrature::wts4, 4}, // QUAD=4
    {GaussQuadrature::pts5, GaussQuadrature::wts5, 5}  // QUAD=5
};

// Global variables of the quadrature and precomputed constants
// Quadrature global variables and constants
const double* Q_POINTS;
const double* Q_WEIGHTS;
int Q_SIZE;
/// Initialize quadrature informations
void init_quadrature_globals() {
    const QuadratureData& qd = quadratureData[QUAD];
    Q_POINTS  = qd.points;
    Q_WEIGHTS = qd.weights;
    Q_SIZE    = qd.quad_size;
}
// Precomputed values
double alpha_p2;
double xmin_pow;
/// Initialize precomputed constants
void init_precomputed_constants() {
    alpha_p2 = ALPHA + 2.0;
    xmin_pow = std::pow(X_MIN, alpha_p2);
}

/**
 * @brief Computes Legendre polynomial P_i(x) for i <= 4 (unrolled version)
 * @param i Polynomial order (0–4)
 * @param x Evaluation point
 * @return P_i(x)
 */
inline double legendre(int i, double x) noexcept {
    switch (i) {
        case 0: return 1.0;
        case 1: return x;
        case 2: return 0.5 * (3.0 * x * x - 1.0);
        case 3: return 0.5 * (5.0 * x * x * x - 3.0 * x);
        case 4: return (35.0 * x * x * x * x - 30.0 * x * x + 3.0) / 8.0;
        default: return 0.0;  // Safe fallback
    }
}

/**
 * @brief Computes derivative of Legendre polynomial P_i'(x) for i <= 4 (unrolled version)
 * @param i Polynomial order (0–4)
 * @param x Evaluation point
 * @return P_i'(x)
 */
inline double Dlegendre(int i, double x) noexcept {
    switch (i) {
        case 0: return 0.0;
        case 1: return 1.0;
        case 2: return 3.0 * x;
        case 3: return 0.5 * (15.0 * x * x - 3.0);
        case 4: return 0.5 * (35.0 * x * x * x - 15.0 * x);
        default: return 0.0;  // Safe fallback
    }
}

/**
 * @brief Initialize precomputed Legendre polynomial table laqp[alpha][i]
 * @return 2D vector of size quad_size x (POLY_ORDER + 1)
 */
std::vector<std::vector<double> > init_laqp() {
    const double* points  = Q_POINTS;
    const double* weights = Q_WEIGHTS;
    const int quad_size   = Q_SIZE;
   
    std::vector<std::vector<double> > laqp(quad_size, std::vector<double>(POLY_ORDER + 1));

    for (int alpha = 0; alpha < quad_size; ++alpha)
        for (int i = 0; i <= POLY_ORDER; ++i)
            laqp[alpha][i] = legendre(i, points[alpha]);

    return laqp;
}

/**
 *@brief Initialize precomputed Legendre polynomial table Dlaqp[alpha][i]
 * @return 2D vector of size quad_size x (POLY_ORDER + 1)
 */
std::vector<std::vector<double> > init_Dlaqp() {
    const double* points  = Q_POINTS;
    const double* weights = Q_WEIGHTS;
    const int quad_size   = Q_SIZE;
    
    std::vector<std::vector<double> > Dlaqp(quad_size, std::vector<double>(POLY_ORDER + 1));

    for (int alpha = 0; alpha < quad_size; ++alpha)
        for (int i = 0; i <= POLY_ORDER; ++i)
            Dlaqp[alpha][i] = Dlegendre(i, points[alpha]);

    return Dlaqp;
}

/**
 * @brief Initial mass density function
 * @param x Mass coordinate
 * @return Initial density value
 */
double g_0(double x) {
    return x * std::exp(-x); // Example: g_0(x) = x * exp(-x)
}

/**
 * @brief Root mean square relative velocity function
 * @param u First mass
 * @param v Second mass
 * @return RMS relative velocity
 */
double rms_v(double u, double v) {
    return 0.5; // Constant velocity for simplicity
}

/**
 * @brief Coagulation kernel function
 * @param u First mass
 * @param v Second mass
 * @param coagulationKernel Selected kernel type
 * @return Kernel value
 */
double coagkernel(double u, double v, CoagulationKernel coagulationKernel) {
    switch (coagulationKernel) {
        case COAG_CONSTANT:
            return 1.0;
        case COAG_ADDITIVE:
            return u + v;
        case COAG_MULTIPLICATIVE:
            return u * v;
        case COAG_BALLISTIC: {
            double cu = std::cbrt(u);
            double cv = std::cbrt(v);
            double sum = cu + cv;
            return M_PI * sum * sum;
        }
        case COAG_DETER:{
            double velo = rms_v(u, v);
            if (velo >= VELO_B) {
                return 0;
            } else {
                double cu = std::cbrt(u);
                double cv = std::cbrt(v);
                double sum = cu + cv;
                return M_PI * sum * sum * velo;
            }
        }
        case COAG_PROB:{
            double velo = rms_v(u, v);
            double term1 = std::sqrt(3.0 / (2.0 * velo * velo));
            double erfResult = std::erf(term1 * VELO_B);

            double term2 = VELO_B * std::sqrt(6.0 / (M_PI * velo * velo));
            double expArg = -3.0 * VELO_B * VELO_B / (2.0 * velo * velo);
            double expResult = std::exp(expArg);
            double cu = std::cbrt(u);
            double cv = std::cbrt(v);
            double sum = cu + cv;
            double sigma_v = M_PI * sum * sum * velo;
            return (erfResult - term2 * expResult) * sigma_v;
        }
        default:
            throw std::invalid_argument("Unknown coagulation kernel type");
    }
}

/**
 * @brief Fragmentation kernel function
 * @param u First mass
 * @param v Second mass
 * @param fragmentationKernel Selected kernel type
 * @return Kernel value
 */
double fragkernel(double u, double v, FragmentationKernel fragmentationKernel) {
    switch (fragmentationKernel) {
        case FRAG_CONSTANT:
            return 1.0;
        case FRAG_ADDITIVE:
            return u + v;
        case FRAG_MULTIPLICATIVE:
            return u * v;
        case FRAG_BALLISTIC:{
            double cu = std::cbrt(u);
            double cv = std::cbrt(v);
            double sum = cu + cv;
            return M_PI * sum * sum;
        }
        case FRAG_DETER:{
            double velo = rms_v(u, v);
            if (velo >= VELO_F) {
                double cu = std::cbrt(u);
                double cv = std::cbrt(v);
                double sum = cu + cv;
                return M_PI * sum * sum * velo;
            } else {
                return 0;
            }
        }
        case FRAG_PROB:{
            double velo = rms_v(u, v);
            double term1 = std::sqrt(3.0 / (2.0 * velo * velo));
            double erfResult = std::erf(term1 * VELO_F);

            double term2 = VELO_F * std::sqrt(6.0 / (M_PI * velo * velo));
            double expArg = -3.0 * VELO_F * VELO_F / (2.0 * velo * velo);
            double expResult = std::exp(expArg);
            double cu = std::cbrt(u);
            double cv = std::cbrt(v);
            double sum = cu + cv;
            double sigma_v = M_PI * sum * sum * velo;
            return (1.0 - erfResult + term2 * expResult) * sigma_v;
        }
        default:
            throw std::invalid_argument("Unknown fragmentation kernel type");
    }
}

/**
 * @brief Calculate total fragment mass after collision
 * @param u First mass
 * @param v Second mass
 * @return Fragment mass
 */
double m_frag(double u, double v) {
    // REFs: Eqs. 34 in Li+2022;
    // Modified with two smooth cosine transitions to prevent numerical oscillations.
    
    double ratio = std::max(u, v) / std::min(u, v);
    double total_mass = u + v;
    
    if (ratio <= 10.0) {
        // Fully destructive
        return total_mass;
    } 
    
    // For all regimes where ratio > 10, calculate the common base factor once
    double base_factor = total_mass / (1.0 + ratio);
    
    if (ratio <= 12.0) {
        // Transition: Destructive -> Cratering
        // f1(q) = 6.5 + 4.5 * cos((q - 10) / 2 * pi)
        double g = 6.5 + 4.5 * std::cos((ratio - 10.0) / 2.0 * M_PI);
        return base_factor * g;
        
    } else if (ratio <= 15.0) {
        // Cratering
        return base_factor * 2.0;
        
    } else if (ratio <= 50.0) {
        // Transition: Cratering -> Mass Transfer
        // f(q) = 1.45 + 0.55 * cos((q - 15) / 35 * pi)
        double f = 1.45 + 0.55 * std::cos((ratio - 15.0) / 35.0 * M_PI);
        return base_factor * f;
        
    } else {
        // Mass transfer
        return base_factor * 0.9;
    }
}

/**
 * @brief Helper function for fragmentation integral
 * @param u First mass
 * @param v Second mass
 * @param x Mass coordinate
 * @return Integral value
 */
double inte1(double u, double v, double x) {
    double integral = 0.0;
    double mfrag = m_frag(u, v);
    double mleft = u + v - mfrag;
    if (x < mleft) {
        integral += mleft;
    }
    double lower = std::min(mfrag, x);
    integral += mfrag * (std::pow(mfrag, alpha_p2) - std::pow(lower, alpha_p2)) / (std::pow(mfrag, alpha_p2) - xmin_pow);
    return integral;
}

/**
 * @brief Helper function for fragmentation integral
 * @param u First mass
 * @param v Second mass
 * @param x Mass coordinate
 * @return Integral value
 */
double inte2(double u, double v, double x) {
    double integral = 0.0;
    double mfrag = m_frag(u, v);
    double mleft = u + v - mfrag;
    if (x >= mleft) {
        integral += mleft;
    }
    double upper = std::min(mfrag, x);
    integral += mfrag * (std::pow(upper, alpha_p2) - xmin_pow) / (std::pow(mfrag, alpha_p2) - xmin_pow);
    return integral;
}

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
void initialize_mass_density(std::vector<std::vector<double> >& g, std::vector<double>& h, std::vector<double>& x, std::vector<double>& cell) {
    const QuadratureData& qd = quadratureData[5];
    const double* points = qd.points;
    const double* weights = qd.weights;
    const int quad_size = qd.quad_size;
    for (int j = 0; j < NUM_BINS; ++j) {
        for (int i = 0; i < POLY_ORDER + 1; ++i) {
            double integral = 0.0;
            for (int k = 0; k < quad_size; ++k) {
                integral += weights[k] * g_0(0.5 * h[j] * points[k] + x[j]) * legendre(i, points[k]);
            }
            g[j][i] = integral * (2.0 * i + 1.0) / 2.0;
        }
    }
    // Apply floor value
    for (int j = 0; j < NUM_BINS; ++j) {
        if (std::abs(g[j][0] * h[j]) < 1e-15) {
            g[j][0] = 1e-15 * std::copysign(1.0, g[j][0]) / h[j];
            for (int i = 1; i < POLY_ORDER + 1; i++) {
                g[j][i] = 0.0;
            }
        }
    }
}

/**
 * @brief Calculate cell-averaged mass density
 * @param g DG coefficients
 * @param g_bar Output: Cell averages
 */
void cal_cell_average(std::vector<std::vector<double> >& g, std::vector<double>& g_bar) {
    for (int j = 0; j < NUM_BINS; ++j) {
        g_bar[j] = g[j][0]; // Zeroth coefficient is average
        if (g_bar[j] < 0) {
            std::cout << "averaged mass density is nagative at j = " << j << std::endl;
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
void cal_gamma(std::vector<std::vector<double> >& g, std::vector<double>& h, std::vector<double>& x,
               std::vector<double>& cell, std::vector<double>& g_bar, std::vector<double>& gamma) {
    cal_cell_average(g, g_bar);
    if (POLY_ORDER == 0) {
        std::fill(gamma.begin(), gamma.end(), 1.0);
        return;
    }
    for (int j = 0; j < NUM_BINS; ++j) {
        std::vector<double> extrema;
        // Boundary points
        extrema.push_back(-1.0);
        extrema.push_back(1.0);

        // Find extrema in the cell
        if (POLY_ORDER >= 2) {
            if (g[j][2] != 0) {
                double x_value = -g[j][1] / (3.0 * g[j][2]);
                if (x_value > -1.0 && x_value < 1.0) {
                    extrema.push_back(x_value);
                }
            }
        }
        if (POLY_ORDER >= 3) {
            if (g[j][3] != 0) {
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
        }
        if (POLY_ORDER >= 4) {
            if (g[j][4] != 0) {
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
        }
        double g_min = std::numeric_limits<double>::max();
        for (int k = 0; k < extrema.size(); ++k) {
            double x_val = extrema[k];
            double poly_val = 0.0;
            for (int i = 0; i <= POLY_ORDER; ++i) {
                poly_val += g[j][i] * legendre(i, x_val);
            }
            if (poly_val < g_min) {
                g_min = poly_val;
            }
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
void cal_p(std::vector<std::vector<double> >& g, std::vector<double>& gamma, std::vector<std::vector<double> >& p) {
    for (int j = 0; j < NUM_BINS; ++j) {
        p[j][0] = g[j][0]; 
        for (int i = 1; i < POLY_ORDER + 1; ++i) {
            p[j][i] = g[j][i] * gamma[j];
        }
    }
}

/**
 * @brief Finds the bin index for a given Value
 * @param cell Bin boundaries
 * @param val Value to locate
 * @return Bin index
 */
inline int find_bin_fast(const std::vector<double>& cell, double val) noexcept {
    auto it = std::lower_bound(cell.begin(), cell.end(), val);
    return static_cast<int>(it - cell.begin()) - 1;
}

/**
 * @brief Calculates coagulation flux at cell boundaries (non-conservative version)
 * @param p Scaled polynomial coefficients
 * @param h Cell widths
 * @param x Cell centers
 * @param cell Bin boundaries
 * @param j Boundary index
 * @return Flux value at boundary xjp1o2
 */
double cal_flux_coag_nc(std::vector<std::vector<double> >& p, std::vector<double>& h, std::vector<double>& x, std::vector<double>& cell, int j, const std::vector<std::vector<double> >& laqp) {
    if (j == -1) return 0.0;

    const double* points  = Q_POINTS;
    const double* weights = Q_WEIGHTS;
    const int quad_size   = Q_SIZE;
    
    double flux = 0.0;
    for (int l = 0; l < j + 1; l++) {
        double flux_l = 0.0;
        for (int alpha = 0; alpha < quad_size; alpha++) {
            double g_l_alpha = 0.0;
            double Gamma_l_alpha = 0.0;
            for (int i = 0; i < POLY_ORDER + 1; i++) {
                g_l_alpha += p[l][i] * laqp[alpha][i];
            }
            double x_l_alpha = x[l] + 0.5 * h[l] * points[alpha];
            double aj = cell[j + 1] - x_l_alpha + cell[0];
            int A = find_bin_fast(cell, aj);
            for (int m = A; m < NUM_BINS; m++) {
                double am = cell[m];
                double bm = cell[m + 1];
                if (m == A) {
                    am = aj;
                }
                double flux_m = 0.0;
                double halfplusm = 0.5 * (bm + am);
                double halfminusm = 0.5 * (bm - am);
                double invhm = 2.0 / h[m];
                for (int beta = 0; beta < quad_size; beta++) {
                    double x_m_beta = halfplusm + halfminusm * points[beta];
                    double xi_m_beta = (x_m_beta - x[m]) * invhm;
                    double g_m_beta = 0.0;
                    for (int i = 0; i < POLY_ORDER + 1; i++) {
                        g_m_beta += p[m][i] * legendre(i, xi_m_beta);
                    }
                    flux_m += weights[beta] * coagkernel(x_l_alpha, x_m_beta, coagulationKernel) * g_m_beta / x_m_beta;
                }
                Gamma_l_alpha += halfminusm * flux_m;
            }
            flux_l += weights[alpha] * g_l_alpha * Gamma_l_alpha;
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
double cal_flux_integral_coag_nc(std::vector<std::vector<double> >& p, std::vector<double>& h, std::vector<double>& x, std::vector<double>& cell, int j, int i, const std::vector<std::vector<double> >& Dlaqp) {
    if (i == 0) return 0.0;
    
    const double* points  = Q_POINTS;
    const double* weights = Q_WEIGHTS;
    const int quad_size   = Q_SIZE;

    double flux_integral = 0.0;
    for (int theta = 0;  theta < quad_size; theta++) {
        double flux_j_theta = 0.0;
        double x_j_theta = x[j] + 0.5 * h[j] * points[theta];
        for (int l = 0; l < j + 1; l++) {
            double al = cell[l];
            double bl = cell[l + 1];
            if (l == j) {
                bl = x_j_theta;
            }
            double halfplusl = 0.5 * (bl + al);
            double halfminusl = 0.5 * (bl - al);
            double invhl = 2.0 / h[l];
            double flux_l = 0.0;
            for (int alpha = 0; alpha < quad_size; alpha++) {
                double x_l_alpha = halfplusl + halfminusl * points[alpha];
                double xi_l_alpha = (x_l_alpha - x[l]) * invhl;
                double g_l_alpha = 0.0;
                double Gamma_l_alpha = 0.0;
                for (int k = 0; k < POLY_ORDER + 1; k++) {
                    g_l_alpha += p[l][k] * legendre(k, xi_l_alpha);
                }
                double aj = x_j_theta - x_l_alpha + cell[0];
                int A = find_bin_fast(cell, aj);
                for (int m = A; m < NUM_BINS; m++) {
                    double am = cell[m];
                    double bm = cell[m + 1];
                    if (m == A) {
                        am = aj;
                    }
                    double flux_m = 0.0;
                    double halfplusm = 0.5 * (bm + am);
                    double halfminusm = 0.5 * (bm - am);
                    double invhm = 2.0 / h[m];
                    for (int beta = 0; beta < quad_size; beta++) {
                        double x_m_beta = halfplusm + halfminusm * points[beta];
                        double xi_m_beta = (x_m_beta - x[m]) * invhm;
                        double g_m_beta = 0.0;
                        for (int k = 0; k < POLY_ORDER + 1; k++) {
                            g_m_beta += p[m][k] * legendre(k, xi_m_beta);
                        }
                        flux_m += weights[beta] * coagkernel(x_l_alpha, x_m_beta, coagulationKernel) * g_m_beta / x_m_beta;
                    }
                    Gamma_l_alpha += halfminusm * flux_m;
                }
                flux_l += weights[alpha] * g_l_alpha * Gamma_l_alpha;
            }
            flux_j_theta += halfminusl * flux_l;
        }
        flux_integral += weights[theta] * Dlaqp[theta][i] * flux_j_theta;
    }
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

    const double* points  = Q_POINTS;
    const double* weights = Q_WEIGHTS;
    const int quad_size   = Q_SIZE;

    double flux = 0.0;
    for (int l = 0; l < j + 1; l++) {
        double flux_l = 0.0;
        for (int alpha = 0; alpha < quad_size; alpha++) {
            double x_l_alpha = x[l] + 0.5 * h[l] * points[alpha];
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
                for (int beta = 0; beta < quad_size; beta++) {
                    double x_m_beta = halfplusm + halfminusm * points[beta];
                    double xi_m_beta = (x_m_beta - x[m]) * invhm;
                    double g_m_beta = 0.0;
                    for (int i = 0; i < POLY_ORDER + 1; i++) {
                        g_m_beta += p[m][i] * legendre(i, xi_m_beta);
                    }
                    flux_m += weights[beta] * coagkernel(x_l_alpha, x_m_beta, coagulationKernel) * g_m_beta / x_m_beta;
                }
                Gamma_l_alpha += halfminusm * flux_m;
            }
            flux_l += weights[alpha] * g_l_alpha * Gamma_l_alpha;
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

    const double* points  = Q_POINTS;
    const double* weights = Q_WEIGHTS;
    const int quad_size   = Q_SIZE;

    double flux_integral = 0.0;
    for (int theta = 0;  theta < quad_size; theta++) {
        double flux_j_theta = 0.0;
        double x_j_theta = x[j] + 0.5 * h[j] * points[theta];
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
            for (int alpha = 0; alpha < quad_size; alpha++) {
                double x_l_alpha = halfplusl + halfminusl * points[alpha];
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
                    for (int beta = 0; beta < quad_size; beta++) {
                        double x_m_beta = halfplusm + halfminusm * points[beta];
                        double xi_m_beta = (x_m_beta - x[m]) * invhm;
                        double g_m_beta = 0.0;
                        for (int k = 0; k < POLY_ORDER + 1; k++) {
                            g_m_beta += p[m][k] * legendre(k, xi_m_beta);
                        }
                        flux_m += weights[beta] * coagkernel(x_l_alpha, x_m_beta, coagulationKernel) * g_m_beta / x_m_beta;
                    }
                    Gamma_l_alpha += halfminusm * flux_m;
                }
                flux_l += weights[alpha] * g_l_alpha * Gamma_l_alpha;
            }
            flux_j_theta += halfminusl * flux_l;
        }
        flux_integral += weights[theta] * Dlaqp[theta][i] * flux_j_theta;
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

    const double* points  = Q_POINTS;
    const double* weights = Q_WEIGHTS;
    const int quad_size   = Q_SIZE;

    double flux1 = 0.0;
    for (int l = 0; l < j + 1; l++) {
        double flux_l = 0.0;
        for (int alpha = 0; alpha < quad_size; alpha++) {
            double x_l_alpha = x[l] + 0.5 * h[l] * points[alpha];
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
                for (int beta = 0; beta < quad_size; beta++) {
                    double x_m_beta = halfplusm + halfminusm * points[beta];
                    double xi_m_beta = (x_m_beta - x[m]) * invhm;
                    double g_m_beta = 0.0;
                    for (int i = 0; i < POLY_ORDER + 1; i++) {
                        g_m_beta += p[m][i] * legendre(i, xi_m_beta);
                    }
                    double integral = inte1(x_l_alpha, x_m_beta, cell[j + 1]);
                    flux_m += weights[beta] * g_m_beta / x_m_beta / (x_m_beta + x_l_alpha) * integral * fragkernel(x_l_alpha, x_m_beta, fragmentationKernel);
                }
                Gamma_l_alpha += halfminusm * flux_m;
            }
            flux_l += weights[alpha] * g_l_alpha * Gamma_l_alpha;
        }
        flux1 += h[l] * 0.5 * flux_l;
    }

    double flux2 = 0.0;
    for (int l = j + 1; l < NUM_BINS; l++) {
        double flux_l = 0.0;
        for (int alpha = 0; alpha < quad_size; alpha++) {
            double x_l_alpha = x[l] + 0.5 * h[l] * points[alpha];
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
                for (int beta = 0; beta < quad_size; beta++) {
                    double x_m_beta = halfplusm + halfminusm * points[beta];
                    double xi_m_beta = (x_m_beta - x[m]) * invhm;
                    double g_m_beta = 0.0;
                    for (int i = 0; i < POLY_ORDER + 1; i++) {
                        g_m_beta += p[m][i] * legendre(i, xi_m_beta);
                    }
                    double integral = inte2(x_l_alpha, x_m_beta, cell[j + 1]);
                    flux_m += weights[beta] * g_m_beta / x_m_beta / (x_m_beta + x_l_alpha) * integral * fragkernel(x_l_alpha, x_m_beta, fragmentationKernel);
                }
                Gamma_l_alpha += halfminusm * flux_m;
            }
            flux_l += weights[alpha] * g_l_alpha * Gamma_l_alpha;
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

    const double* points  = Q_POINTS;
    const double* weights = Q_WEIGHTS;
    const int quad_size   = Q_SIZE;

    double flux1 = 0.0;
    for (int l = 0; l < NUM_BINS; l++) {
        double flux_l = 0.0;
        for (int alpha = 0; alpha < quad_size; alpha++) {
            double x_l_alpha = x[l] + 0.5 * h[l] * points[alpha];
            double g_l_alpha = 0.0;
            for (int k = 0; k < POLY_ORDER + 1; k++) {
                g_l_alpha += p[l][k] * laqp[alpha][i];
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
                for (int beta = 0; beta < quad_size; beta++) {
                    double x_m_beta = halfplusm + halfminusm * points[beta];
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
                        for (int theta = 0;  theta < quad_size; theta++) {
                            double xi_j_theta = halfplusl + halfminusl * points[theta];
                            double x_j_theta = 0.5 * h[j] * xi_j_theta + x[j];
                            inner2 += halfminusl * weights[theta] * Dlegendre(i, xi_j_theta) * std::pow(x_j_theta, alpha_p2);
                        }
                        integral += mfrag / (std::pow(mfrag, alpha_p2) - xmin_pow) * (inner1 - inner2);
                    }
                    flux_m += weights[beta] * g_m_beta / x_m_beta / (x_m_beta + x_l_alpha) * fragkernel(x_l_alpha, x_m_beta, fragmentationKernel) * integral;
                }
                Gamma_l_alpha += halfminusm * flux_m;
            }
            flux_l += weights[alpha] * g_l_alpha * Gamma_l_alpha;
        }
        flux1 += h[l] * 0.5 * flux_l;
    }

    double flux2 = 0.0;
    for (int l = 0; l < NUM_BINS; l++) {
        double flux_l = 0.0;
        for (int alpha = 0; alpha < quad_size; alpha++) {
            double x_l_alpha = x[l] + 0.5 * h[l] * points[alpha];
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
                for (int beta = 0; beta < quad_size; beta++) {
                    double x_m_beta = halfplusm + halfminusm * points[beta];
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
                        for (int theta = 0;  theta < quad_size; theta++) {
                            double xi_j_theta = halfplusl + halfminusl * points[theta];
                            double x_j_theta = 0.5 * h[j] * xi_j_theta + x[j];
                            inner2 += halfminusl * weights[theta] * Dlegendre(i, xi_j_theta) * std::pow(x_j_theta, alpha_p2);
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
                        for (int theta = 0;  theta < quad_size; theta++) {
                            double xi_j_theta = halfplusl + halfminusl * points[theta];
                            double x_j_theta = 0.5 * h[j] * xi_j_theta + x[j];
                            inner2 += halfminusl * weights[theta] * Dlegendre(i, xi_j_theta) * std::pow(x_j_theta, alpha_p2);
                        }
                        integral += mfrag / (std::pow(mfrag, alpha_p2) - xmin_pow) * (inner1 - inner2);
                    }
                    flux_m += weights[beta] * g_m_beta / x_m_beta / (x_m_beta + x_l_alpha) * fragkernel(x_l_alpha, x_m_beta, fragmentationKernel) * integral;
                }
                Gamma_l_alpha += halfminusm * flux_m;
            }
            flux_l += weights[alpha] * g_l_alpha * Gamma_l_alpha;
        }
        flux2 += h[l] * 0.5 * flux_l;
    }

    double flux3 = 0.0;
    for (int l = 0; l < NUM_BINS; l++) {
        double flux_l = 0.0;
        for (int alpha = 0; alpha < quad_size; alpha++) {
            double x_l_alpha = x[l] + 0.5 * h[l] * points[alpha];
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
                for (int beta = 0; beta < quad_size; beta++) {
                    double x_m_beta = halfplusm + halfminusm * points[beta];
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
                    flux_m += weights[beta] * g_m_beta * mleft / x_m_beta / (x_m_beta + x_l_alpha) * fragkernel(x_l_alpha, x_m_beta, fragmentationKernel) * integral;
                }
                Gamma_l_alpha += halfminusm * flux_m;
            }
            flux_l += weights[alpha] * g_l_alpha * Gamma_l_alpha;
        }
        flux3 += h[l] * 0.5 * flux_l;
    }

    double flux4 = 0.0;
    for (int l = 0; l < NUM_BINS; l++) {
        double flux_l = 0.0;
        for (int alpha = 0; alpha < quad_size; alpha++) {
            double x_l_alpha = x[l] + 0.5 * h[l] * points[alpha];
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
                for (int beta = 0; beta < quad_size; beta++) {
                    double x_m_beta = halfplusm + halfminusm * points[beta];
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
                    flux_m += weights[beta] * g_m_beta * mleft / x_m_beta / (x_m_beta + x_l_alpha) * fragkernel(x_l_alpha, x_m_beta, fragmentationKernel) * integral;
                }
                Gamma_l_alpha += halfminusm * flux_m;
            }
            flux_l += weights[alpha] * g_l_alpha * Gamma_l_alpha;
        }
        flux4 += h[l] * 0.5 * flux_l;
    }
    double flux_integral = flux1 - flux2 + flux3 - flux4;
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

    const double* points  = Q_POINTS;
    const double* weights = Q_WEIGHTS;
    const int quad_size   = Q_SIZE;

    double flux1 = 0.0;
    for (int l = 0; l < j + 1; l++) {
        double flux_l = 0.0;
        for (int alpha = 0; alpha < quad_size; alpha++) {
            double x_l_alpha = x[l] + 0.5 * h[l] * points[alpha];
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
                for (int beta = 0; beta < quad_size; beta++) {
                    double x_m_beta = halfplusm + halfminusm * points[beta];
                    double xi_m_beta = (x_m_beta - x[m]) * invhm;
                    double g_m_beta = 0.0;
                    for (int i = 0; i < POLY_ORDER + 1; i++) {
                        g_m_beta += p[m][i] * legendre(i, xi_m_beta);
                    }
                    double mfrag = m_frag(x_l_alpha, x_m_beta);
                    double integral = inte1(x_l_alpha, x_m_beta, cell[j + 1]);
                    flux_m += weights[beta] * g_m_beta / x_m_beta / (x_m_beta + x_l_alpha) * integral * fragkernel(x_l_alpha, x_m_beta,fragmentationKernel);
                }
                Gamma_l_alpha += halfminusm * flux_m;
            }
            flux_l += weights[alpha] * g_l_alpha * Gamma_l_alpha;
        }
        flux1 += h[l] * 0.5 * flux_l;
    }

    double flux2 = 0.0;
    for (int l = j + 1; l < NUM_BINS; l++) {
        double flux_l = 0.0;
        for (int alpha = 0; alpha < quad_size; alpha++) {
            double x_l_alpha = x[l] + 0.5 * h[l] * points[alpha];
            double g_l_alpha = 0.0;
            for (int i = 0; i < POLY_ORDER + 1; i++) {
                g_l_alpha += p[l][i] * laqp[alpha][i];
            }
            double Gamma_l_alpha = 0.0;
            for (int m = 0; m < NUM_BINS; m++) {
                double flux_m = 0.0;
                for (int beta = 0; beta < quad_size; beta++) {
                    double x_m_beta = x[m] + 0.5 * h[m] * points[beta];
                    double g_m_beta = 0.0;
                    for (int i = 0; i < POLY_ORDER + 1; i++) {
                        g_m_beta += p[m][i] * legendre(i, points[beta]);
                    }
                    double mfrag = m_frag(x_l_alpha, x_m_beta);
                    double integral = inte2(x_l_alpha, x_m_beta, cell[j + 1]);
                    flux_m += weights[beta] * g_m_beta / x_m_beta / (x_m_beta + x_l_alpha) * integral * fragkernel(x_l_alpha, x_m_beta, fragmentationKernel);
                }
                Gamma_l_alpha += 0.5 * h[m] * flux_m;
            }
            flux_l += weights[alpha] * g_l_alpha * Gamma_l_alpha;
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

    const double* points  = Q_POINTS;
    const double* weights = Q_WEIGHTS;
    const int quad_size   = Q_SIZE;

    double flux1 = 0.0;
    for (int l = 0; l < NUM_BINS; l++) {
        double flux_l = 0.0;
        for (int alpha = 0; alpha < quad_size; alpha++) {
            double x_l_alpha = x[l] + 0.5 * h[l] * points[alpha];
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
                for (int beta = 0; beta < quad_size; beta++) {
                    double x_m_beta = halfplusm + halfminusm * points[beta];
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
                        for (int theta = 0;  theta < quad_size; theta++) {
                            double xi_j_theta = halfplusl + halfminusl * points[theta];
                            double x_j_theta = 0.5 * h[j] * xi_j_theta + x[j];
                            inner2 += halfminusl * weights[theta] * Dlegendre(i, xi_j_theta) * std::pow(x_j_theta, alpha_p2);
                        }
                        integral += mfrag / (std::pow(mfrag, alpha_p2) - xmin_pow) * (inner1 - inner2);
                    }
                    flux_m += weights[beta] * g_m_beta / x_m_beta / (x_m_beta + x_l_alpha) * fragkernel(x_l_alpha, x_m_beta, fragmentationKernel) * integral;
                }
                Gamma_l_alpha += halfminusm * flux_m;
            }
            flux_l += weights[alpha] * g_l_alpha * Gamma_l_alpha;
        }
        flux1 += h[l] * 0.5 * flux_l;
    }

    double flux2 = 0.0;
    for (int l = 0; l < NUM_BINS; l++) {
        double flux_l = 0.0;
        for (int alpha = 0; alpha < quad_size; alpha++) {
            double x_l_alpha = x[l] + 0.5 * h[l] * points[alpha];
            double g_l_alpha = 0.0;
            for (int k = 0; k < POLY_ORDER + 1; k++) {
                g_l_alpha += p[l][k] * laqp[alpha][k];
            }
            double Gamma_l_alpha = 0.0;
            for (int m = 0; m < NUM_BINS; m++) {
                double flux_m = 0.0;
                for (int beta = 0; beta < quad_size; beta++) {
                    double x_m_beta = x[m] + 0.5 * h[m] * points[beta];
                    double g_m_beta = 0.0;
                    for (int k = 0; k < POLY_ORDER + 1; k++) {
                        g_m_beta += p[m][k] * laqp[beta][k];
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
                        for (int theta = 0;  theta < quad_size; theta++) {
                            double xi_j_theta = halfplusl + halfminusl * points[theta];
                            double x_j_theta = 0.5 * h[j] * xi_j_theta + x[j];
                            inner2 += halfminusl * weights[theta] * Dlegendre(i, xi_j_theta) * std::pow(x_j_theta, alpha_p2);
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
                        for (int theta = 0;  theta < quad_size; theta++) {
                            double xi_j_theta = halfplusl + halfminusl * points[theta];
                            double x_j_theta = 0.5 * h[j] * xi_j_theta + x[j];
                            inner2 += halfminusl * weights[theta] * Dlegendre(i, xi_j_theta) * std::pow(x_j_theta, alpha_p2);
                        }
                        integral += mfrag / (std::pow(mfrag, alpha_p2) - xmin_pow) * (inner1 - inner2);
                    }
                    flux_m += weights[beta] * g_m_beta / x_m_beta / (x_m_beta + x_l_alpha) * fragkernel(x_l_alpha, x_m_beta, fragmentationKernel) * integral;
                }
                Gamma_l_alpha += 0.5 * h[m] * flux_m;
            }
            flux_l += weights[alpha] * g_l_alpha * Gamma_l_alpha;
        }
        flux2 += h[l] * 0.5 * flux_l;
    }

    double flux3 = 0.0;
    for (int l = 0; l < NUM_BINS; l++) {
        double flux_l = 0.0;
        for (int alpha = 0; alpha < quad_size; alpha++) {
            double x_l_alpha = x[l] + 0.5 * h[l] * points[alpha];
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
                for (int beta = 0; beta < quad_size; beta++) {
                    double x_m_beta = halfplusm + halfminusm * points[beta];
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
                    flux_m += weights[beta] * g_m_beta * mleft / x_m_beta / (x_m_beta + x_l_alpha) * fragkernel(x_l_alpha, x_m_beta, fragmentationKernel) * integral;
                }
                Gamma_l_alpha += halfminusm * flux_m;
            }
            flux_l += weights[alpha] * g_l_alpha * Gamma_l_alpha;
        }
        flux3 += h[l] * 0.5 * flux_l;
    }

    double flux4 = 0.0;
    for (int l = 0; l < NUM_BINS; l++) {
        double flux_l = 0.0;
        for (int alpha = 0; alpha < quad_size; alpha++) {
            double x_l_alpha = x[l] + 0.5 * h[l] * points[alpha];
            double g_l_alpha = 0.0;
            for (int k = 0; k < POLY_ORDER + 1; k++) {
                g_l_alpha += p[l][k] * laqp[alpha][k];
            }
            double Gamma_l_alpha = 0.0;
            for (int m = 0; m < NUM_BINS; m++) {
                double flux_m = 0.0;
                for (int beta = 0; beta < quad_size; beta++) {
                    double x_m_beta = x[m] + 0.5 * h[m] * points[beta];
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
                    flux_m += weights[beta] * g_m_beta * mleft / x_m_beta / (x_m_beta + x_l_alpha) * fragkernel(x_l_alpha, x_m_beta, fragmentationKernel) * integral;
                }
                Gamma_l_alpha += 0.5 * h[m] * flux_m;
            }
            flux_l += weights[alpha] * g_l_alpha * Gamma_l_alpha;
        }
        flux4 += h[l] * 0.5 * flux_l;
    }
    double flux_integral = flux1 - flux2 + flux3 - flux4;
    return flux_integral;
}

/**
 * @brief Calculates adaptive time step based on flux conditions
 * @param p Scaled polynomial coefficients
 * @param h Cell widths
 * @param x Cell centers
 * @param cell Bin boundaries
 * @param flux Flux values (output)
 * @return Calculated time step
 */
double cal_time_step(std::vector<std::vector<double> >& p, std::vector<double>& h, std::vector<double>& x, std::vector<double>& cell, std::vector<double>& flux, const std::vector<std::vector<double>>& laqp) {
    std::vector<double> delta_t(NUM_BINS);
    double flux_diff;
    for (int j = 0; j < NUM_BINS; j++) {
        if (flux_mode == "coag_frag") {
            flux[j] = cal_flux_coag_nc(p, h, x, cell, j-1, laqp)
                    + cal_flux_frag_nc(p, h, x, cell, j-1, laqp);
        }
        else if (flux_mode == "pure_frag") {
            flux[j] = cal_flux_frag_nc(p, h, x, cell, j-1, laqp);
        }
        else if (flux_mode == "pure_coag") {
            flux[j] = cal_flux_coag_nc(p, h, x, cell, j-1, laqp);
        }
        if (j == NUM_BINS-1) {
            if (flux_mode == "coag_frag") {
                flux[j+1] = cal_flux_coag_nc(p, h, x, cell, j, laqp)
                          + cal_flux_frag_nc(p, h, x, cell, j, laqp);
            }
            else if (flux_mode == "pure_frag") {
                flux[j+1] = cal_flux_frag_nc(p, h, x, cell, j, laqp);
            }
            else if (flux_mode == "pure_coag") {
                flux[j+1] = cal_flux_coag_nc(p, h, x, cell, j, laqp);
            }
        }
    }
    for (int j = 0; j < NUM_BINS; j++) {
        flux_diff = flux[j] - flux[j+1];
        delta_t[j] = h[j] * p[j][0] / std::abs(flux_diff);
        // for debug use, check the time step is limited by what
        // std::cout << "j: " << j << std::endl;
        // std::cout << "hj " << h[j] << std::endl;
        // std::cout << "gj0 " << p[j][0] << std::endl;
        // std::cout << "flux_diff " << flux_diff << std::endl;
        // std::cout << "delta_tj " << delta_t[j] << std::endl;
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
bool evolve(std::vector<std::vector<double> >& g, std::vector<std::vector<double> >& p, std::vector<double>& h, std::vector<double>& x, std::vector<double>& cell, std::vector<double>& flux, double& dt, const std::vector<std::vector<double>>& laqp, const std::vector<std::vector<double>>& Dlaqp) {
    // Benckmark tests use RK3
    std::vector<std::vector<double> > g1; 
    g1.resize(NUM_BINS, std::vector<double>(POLY_ORDER + 1));
    std::vector<std::vector<double> > g2;
    g2.resize(NUM_BINS, std::vector<double>(POLY_ORDER + 1));
    std::vector<std::vector<double> > g3;
    g3.resize(NUM_BINS, std::vector<double>(POLY_ORDER + 1));
    for (int j = 0; j < NUM_BINS; j++) {
        for (int i = 0; i < POLY_ORDER + 1; i++) {
            double update = 0.0;
            if (flux_mode == "coag_frag") {
                update = cal_flux_integral_coag_nc(p, h, x, cell, j, i, Dlaqp)
                       + cal_flux_integral_frag_nc(p, h, x, cell, j, i, laqp);
            }
            else if (flux_mode == "pure_frag") {
                update = cal_flux_integral_frag_nc(p, h, x, cell, j, i, laqp);
            }
            else if (flux_mode == "pure_coag") {
                update = cal_flux_integral_coag_nc(p, h, x, cell, j, i, Dlaqp);
            }
            // Add flux difference terms
            update += -flux[j+1] + flux[j] * std::pow(-1, i);
            g1[j][i] = p[j][i] + dt * (2.0 * i + 1.0) / h[j] * update;
        }
    }
    
    double useless = cal_time_step(g1, h, x, cell, flux, laqp);
    for (int j = 0; j < NUM_BINS; j++) {
        for (int i = 0; i < POLY_ORDER + 1; i++) {
            double update = 0.0;
            if (flux_mode == "coag_frag") {
                update = cal_flux_integral_coag_nc(g1, h, x, cell, j, i, Dlaqp)
                       + cal_flux_integral_frag_nc(g1, h, x, cell, j, i, laqp);
            }
            else if (flux_mode == "pure_frag") {
                update = cal_flux_integral_frag_nc(g1, h, x, cell, j, i, laqp);
            }
            else if (flux_mode == "pure_coag") {
                update = cal_flux_integral_coag_nc(g1, h, x, cell, j, i, Dlaqp);
            }
            update += -flux[j+1] + flux[j] * std::pow(-1, i);
            g2[j][i] = 3.0 / 4.0 * p[j][i] + 1.0 / 4.0 * (g1[j][i] + dt * (2.0 * i + 1.0) / h[j] * update);
        }
    }

    useless = cal_time_step(g2, h, x, cell, flux, laqp);
    for (int j = 0; j < NUM_BINS; j++) {
        for (int i = 0; i < POLY_ORDER + 1; i++) {
            double update = 0.0;
            if (flux_mode == "coag_frag") {
                update = cal_flux_integral_coag_nc(g2, h, x, cell, j, i, Dlaqp)
                       + cal_flux_integral_frag_nc(g2, h, x, cell, j, i, laqp);
            }
            else if (flux_mode == "pure_frag") {
                update = cal_flux_integral_frag_nc(g2, h, x, cell, j, i, laqp);
            }
            else if (flux_mode == "pure_coag") {
                update = cal_flux_integral_coag_nc(g2, h, x, cell, j, i, Dlaqp);
            }
            update += -flux[j+1] + flux[j] * std::pow(-1, i); 
            g[j][i] = 1.0 / 3.0 * p[j][i] + 2.0 / 3.0 * (g2[j][i] + dt * (2.0 * i + 1.0) / h[j] * update);
        }
    }
    
    for (int j = 0; j < NUM_BINS; j++) {
        // Apply floor value
        if (g[j][0] * h[j] < 1e-30 && g[j][0] > 0) {
            g[j][0] = 1e-30 / h[j];
            for (int i = 1; i < POLY_ORDER + 1; i++) {
                g[j][i] = 0.0;
            }
        }
        // Check for negative densities
        if (g[j][0] < 0) {
            std::cerr << "Negative density at j=" << j
                      << ": " << g[j][0] << std::endl;
            dt /= 2.0;
            return false;
        }
    }
    return true;
}

/**
 * @brief Records simulation state to file
 * @param filename Output file name
 * @param g Polynomial coefficients
 * @param cell Bin boundaries
 * @param h Cell widths
 * @param x Cell centers
 * @param step Current time step
 * @param T Current simulation time
 */
void record(const std::string& filename, std::vector<std::vector<double> >& g, std::vector<double>& cell,
            std::vector<double>& h, std::vector<double>& x, int step, double T) {
    std::ofstream fp;
    // If step == 0, overwrite the file; else append to it
    if (step == 0) {
        fp.open(filename, std::ios::trunc);  // Overwrite the file
    } else {
        fp.open(filename, std::ios::app);   // Append to the file
    }

    // If it's the first step (step == 0), save cell, h, x, and g
    if (step == 0) {
        fp << "Step: " << step << std::endl;
        fp << "Time: " << T << std::endl;
        for (int i = 0; i < NUM_BINS + 1; i++) {
            fp << cell[i] << std::endl;
        }
        for (int i = 0; i < NUM_BINS; i++) {
            fp << h[i] << std::endl;
        }
        for (int i = 0; i < NUM_BINS; i++) {
            fp << x[i] << std::endl;
        }
        for (int i = 0; i < NUM_BINS; i++) {
            for (int j = 0; j < POLY_ORDER + 1; j++) {
                fp << g[i][j] << " ";
            }
            fp << std::endl;
        }
    }
    // For other steps, first store step, then store g
    else {
        fp << "Step: " << step << std::endl;
        fp << "Time: " << T << std::endl;
        for (int i = 0; i < NUM_BINS; i++) {
            for (int j = 0; j < POLY_ORDER + 1; j++) {
                fp << g[i][j] << " ";
            }
            fp << std::endl;
        }
    }

    fp.close();
}

/**
 * @brief Reads simulation parameters from input file
 * @param filename Input file name
 */
void readInputFile(const std::string& filename) {
    std::ifstream inputFile(filename);
    if (!inputFile.is_open()) {
        std::cerr << "Error: Could not open input file " << filename << std::endl;
        exit(1);
    }

    std::string line;
    while (std::getline(inputFile, line)) {
        std::string key;
        std::string value;
        std::size_t delimiterPos = line.find('=');
        if (delimiterPos != std::string::npos) {
            key = line.substr(0, delimiterPos);
            value = line.substr(delimiterPos + 1);

            // Trim leading and trailing whitespace from key and value
            key.erase(0, key.find_first_not_of(' '));
            key.erase(key.find_last_not_of(' ') + 1);
            value.erase(0, value.find_first_not_of(' '));
            value.erase(value.find_last_not_of(' ') + 1);

            if (key == "NUM_BINS") {
                NUM_BINS = std::stoi(value);
            } else if (key == "POLY_ORDER") {
                POLY_ORDER = std::stoi(value);
                // QUAD = POLY_ORDER + 1;
            } else if (key == "MAX_TIME") {
                MAX_TIME = std::stod(value);
            } else if (key == "X_MIN") {
                X_MIN = std::stod(value);
            } else if (key == "X_MAX") {
                X_MAX = std::stod(value);
            } else if (key == "CFL") {
                CFL = std::stod(value);
            } else if (key == "COAGULATION_KERNEL") {
                if (value == "constant") {
                    coagulationKernel = COAG_CONSTANT;
                } else if (value == "additive") {
                    coagulationKernel = COAG_ADDITIVE;
                } else if (value == "multiplicative") {
                    coagulationKernel = COAG_MULTIPLICATIVE;
                } else if (value == "ballistic") {
                    coagulationKernel = COAG_BALLISTIC;
                } else if (value == "determined") {
                    coagulationKernel = COAG_DETER;
                } else if (value == "probabilistic") {
                    coagulationKernel = COAG_PROB;
                } else {
                    std::cerr << "Error: Unknown COAGULATION_KERNEL value: " << value << std::endl;
                    exit(1);
                }
            } else if (key == "FRAGMENTATION_KERNEL") {
                if (value == "constant") {
                    fragmentationKernel = FRAG_CONSTANT;
                } else if (value == "additive") {
                    fragmentationKernel = FRAG_ADDITIVE;
                } else if (value == "multiplicative") {
                    fragmentationKernel = FRAG_MULTIPLICATIVE;
                } else if (value == "ballistic") {
                    fragmentationKernel = FRAG_BALLISTIC;
                } else if (value == "determined") {
                    fragmentationKernel = FRAG_DETER;
                } else if (value == "probabilistic") {
                    fragmentationKernel = FRAG_PROB;
                } else {
                    std::cerr << "Error: Unknown COAGULATION_KERNEL value: " << value << std::endl;
                    exit(1);
                }
            } else if (key == "VELO_F") { 
                VELO_F = std::stod(value);
            } else if (key == "VELO_B") {   
                VELO_B = std::stod(value);
            } else if (key == "ALPHA") {
                ALPHA = std::stod(value);
            } else if (key == "QUAD") {
                QUAD = std::stoi(value);
            } else if (key == "flux_mode") {
                flux_mode = value;
                // Validate mode
                if (flux_mode != "coag_frag" &&
                    flux_mode != "pure_frag" &&
                    flux_mode != "pure_coag") {
                    std::cerr << "Warning: Invalid flux_mode '" << value
                              << "'. Using default 'coag_frag'." << std::endl;
                    flux_mode = "coag_frag";
                }
            } else {
                std::cerr << "Warning: Unknown key '" << key << "' in input file." << std::endl;
            }
        }
    }
    inputFile.close();
}


/**
 * @brief Main simulation driver
 * @return Exit status
 */
int main() {
    TimeProfile timer;
    readInputFile("input.txt");
    
    // Print configuration
    std::cout << "NUM_BINS: " << NUM_BINS << std::endl;
    std::cout << "POLY_ORDER: " << POLY_ORDER << std::endl;
    std::cout << "QUAD: " << QUAD << std::endl;
    std::cout << "MAX_TIME: " << MAX_TIME << std::endl;
    std::cout << "X_MIN: " << X_MIN << std::endl;
    std::cout << "X_MAX: " << X_MAX << std::endl;
    std::cout << "CFL: " << CFL << std::endl;
    std::cout << "VELO_F: " << VELO_F << std::endl;
    std::cout << "VELO_B: " << VELO_B << std::endl;
    std::cout << "ALPHA: " << ALPHA << std::endl;
    std::cout << "Flux Mode: " << flux_mode << std::endl;
    // Initialize quad
    init_quadrature_globals();
    // Precomputed constants
    init_precomputed_constants(); 
    std::cout << "alpha_p2 = " << alpha_p2 << ", xmin_pow = " << xmin_pow << std::endl;

    // Initialize laqp and Dlaqp
    std::vector<std::vector<double>> laqp = init_laqp();
    std::vector<std::vector<double>> Dlaqp = init_Dlaqp();

    // Initialize grid
    std::vector<double> h(NUM_BINS), x(NUM_BINS), cell(NUM_BINS + 1);
    initialize_mass_bin(h, x, cell);

    // Initialize density
    std::vector<std::vector<double> > g;
    g.resize(NUM_BINS, std::vector<double>(POLY_ORDER + 1));
    std::vector<double> g_bar(NUM_BINS);
    std::vector<double> gamma(NUM_BINS);
    initialize_mass_density(g, h, x, cell);
    cal_gamma(g, h, x, cell, g_bar, gamma);

    std::vector<std::vector<double> > p;
    p.resize(NUM_BINS, std::vector<double>(POLY_ORDER + 1));
    cal_p(g, gamma, p);

    int step = 0;
    double T = 0;
    std::vector<double> flux(NUM_BINS+1);
    // Record initial state at step 0 (overwrite the file)
    record("output.txt", p, cell, h, x, step, T);
    
    // Main loop for simulation
    while (T < MAX_TIME) {
        auto loop_start = std::chrono::high_resolution_clock::now();
        // Calculate time step
        auto t1 = std::chrono::high_resolution_clock::now();
        double dt = cal_time_step(p, h, x, cell, flux, laqp);
        auto t2 = std::chrono::high_resolution_clock::now();
        timer.total_cal_time_step += std::chrono::duration<double>(t2 - t1).count();
        if (T + dt > MAX_TIME) {
            dt = MAX_TIME - T;
        }
        bool success = false;
        while (!success) {
            // Evolve system
            std::cout << "Time step: " << dt << std::endl;
            auto t3 = std::chrono::high_resolution_clock::now();
            success = evolve(g, p, h, x, cell, flux, dt, laqp, Dlaqp);
            auto t4 = std::chrono::high_resolution_clock::now();
            timer.total_evolve += std::chrono::duration<double>(t4 - t3).count();
        }
       
        // Compute gamma and update p
        auto t5 = std::chrono::high_resolution_clock::now();
        cal_gamma(g, h, x, cell, g_bar, gamma);
        auto t6 = std::chrono::high_resolution_clock::now();
        timer.total_cal_gamma += std::chrono::duration<double>(t6 - t5).count();
        
        auto t7 = std::chrono::high_resolution_clock::now();
        cal_p(g, gamma, p);
        auto t8 = std::chrono::high_resolution_clock::now();
        timer.total_cal_p += std::chrono::duration<double>(t8 - t7).count();
        
        // Update time and step
        T += dt;
        step++;
        
        // Record state
        auto t9 = std::chrono::high_resolution_clock::now();
        record("output.txt", p, cell, h, x, step, T);
        auto t10 = std::chrono::high_resolution_clock::now();
        timer.total_record += std::chrono::duration<double>(t10 - t9).count();

        auto loop_end = std::chrono::high_resolution_clock::now();
        timer.total_loop += std::chrono::duration<double>(loop_end - loop_start).count();
        
        std::cout << "Step: " << step
                  << " Time: " << T
                  << " dt: " << dt << std::endl;
    }
     std::cout << "\n--- Total Step " << step << " ---\n"
                  << "cal_time_step: " << timer.total_cal_time_step << " s\n"
                  << "evolve: " << timer.total_evolve << " s\n"
                  << "cal_gamma: " << timer.total_cal_gamma << " s\n"
                  << "cal_p: " << timer.total_cal_p << " s\n"
                  << "record: " << timer.total_record << " s\n"
                  << "Total loop: " << timer.total_loop << " s\n";
    return 0;
}
