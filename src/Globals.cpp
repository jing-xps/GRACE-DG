#include "Globals.hpp"
#include <cmath>

// This file defines global variables and constants used across the coagulation-fragmentation simulation. 

// Definition of global variables
int NUM_BINS = 0;           ///< Number of size bins
int POLY_ORDER = 0;         ///< Polynomial order for DG method (0 to 4)
int QUAD = 0;               ///< Quadrature order (1 to 5)
double MAX_TIME = 0.0;      ///< Maximum simulation time
double X_MIN = 0.0;         ///< Minimum particle mass
double X_MAX = 0.0;         ///< Maximum particle mass
double CFL = 0.0;           ///< CFL number for time-stepping
double VELO_F = 0.0;        ///< Velocity for fragmentation (if applicable)
double VELO_B = 0.0;        ///< Velocity for bouncing (if applicable)
double ALPHA = 0.0;         ///< Fragment distribution power-law index

CoagulationKernel coagulationKernel = COAG_CONSTANT;
FragmentationKernel fragmentationKernel = FRAG_CONSTANT;
std::string flux_mode = "coag_frag";

const double* Q_POINTS = nullptr;
const double* Q_WEIGHTS = nullptr;
int Q_SIZE = 0;
double alpha_p2 = 0.0;
double xmin_pow = 0.0;

/**
 * @namespace GaussQuadrature
 * @brief Contains Gauss quadrature points and weights for different orders
 */
namespace GaussQuadrature {
    const double pts1[] = {0.0};
    const double wts1[] = {2.0};
    const double pts2[] = {-1.0 / 1.732050807568877, 1.0 / 1.732050807568877};
    const double wts2[] = {1.0, 1.0};
    const double pts3[] = {-0.774596669241483, 0.0, 0.774596669241483};
    const double wts3[] = {5.0 / 9.0, 8.0 / 9.0, 5.0 / 9.0};
    const double sqrt30 = 5.477225575051661;
    const double pts4[] = {
        -std::sqrt((3.0 + 2.0 * sqrt30 / 5.0) / 7.0),
        -std::sqrt((3.0 - 2.0 * sqrt30 / 5.0) / 7.0),
        std::sqrt((3.0 - 2.0 * sqrt30 / 5.0) / 7.0),
        std::sqrt((3.0 + 2.0 * sqrt30 / 5.0) / 7.0)
    };
    const double wts4[] = {
        (18.0 - sqrt30) / 36.0, (18.0 + sqrt30) / 36.0,
        (18.0 + sqrt30) / 36.0, (18.0 - sqrt30) / 36.0
    };
    const double sqrt70 = 8.366600265340756;
    const double pts5[] = {
        -std::sqrt(5.0 + 2.0 * sqrt70 / 7.0) / 3.0,
        -std::sqrt(5.0 - 2.0 * sqrt70 / 7.0) / 3.0,
        0.0,
        std::sqrt(5.0 - 2.0 * sqrt70 / 7.0) / 3.0,
        std::sqrt(5.0 + 2.0 * sqrt70 / 7.0) / 3.0
    };
    const double wts5[] = {
        (322.0 - 13.0 * sqrt70) / 900.0, (322.0 + 13.0 * sqrt70) / 900.0,
        128.0 / 225.0, (322.0 + 13.0 * sqrt70) / 900.0, (322.0 - 13.0 * sqrt70) / 900.0
    };
}

/**
 * @struct QuadratureData
 * @brief Stores quadrature points, weights and size
 */
struct QuadratureData {
    const double* points;
    const double* weights;
    int quad_size;
};

/// Quadrature data for orders 1-5
const QuadratureData quadratureData[] = {
    {nullptr, nullptr, 0}, 
    {GaussQuadrature::pts1, GaussQuadrature::wts1, 1}, 
    {GaussQuadrature::pts2, GaussQuadrature::wts2, 2}, 
    {GaussQuadrature::pts3, GaussQuadrature::wts3, 3}, 
    {GaussQuadrature::pts4, GaussQuadrature::wts4, 4}, 
    {GaussQuadrature::pts5, GaussQuadrature::wts5, 5}  
};

/// Initialize quadrature informations
void init_quadrature_globals() {
    const QuadratureData& qd = quadratureData[QUAD];
    Q_POINTS  = qd.points;
    Q_WEIGHTS = qd.weights;
    Q_SIZE    = qd.quad_size;
}

/// Initialize precomputed constants
void init_precomputed_constants() {
    alpha_p2 = ALPHA + 2.0;
    xmin_pow = std::pow(X_MIN, alpha_p2);
}
