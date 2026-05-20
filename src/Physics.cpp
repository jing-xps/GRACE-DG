#include "Physics.hpp"
#include <cmath>
#include <algorithm>
#include <stdexcept>

// This file implements the physics of coagulation and fragmentation processes
// including various kernels and the mass distribution of fragments.

/**
 * @brief Initial mass density function
 * @param x Mass coordinate
 * @return Initial density value
 */
double g_0(double x) {
    return x * std::exp(-x); // Example: g_0(x) = x * exp(-x)
}

double rms_v(double u, double v) {
    (void)u; (void)v;
    return 0.5; // Constant velocity for simplicity
}

/**
 * @brief Coagulation kernel function
 * @param u First mass
 * @param v Second mass
 * @param coagulationKernel Selected kernel type
 * @return Kernel value
 */
double coagkernel(double u, double v, CoagulationKernel kernel) {
    switch (kernel) {
        case COAG_CONSTANT: return 1.0;
        case COAG_ADDITIVE: return u + v;
        case COAG_MULTIPLICATIVE: return u * v;
        case COAG_BALLISTIC: {
            double sum = std::cbrt(u) + std::cbrt(v);
            return M_PI * sum * sum;
        }
        case COAG_DETER: {
            double velo = rms_v(u, v);
            if (velo >= VELO_B) return 0.0;
            double sum = std::cbrt(u) + std::cbrt(v);
            return M_PI * sum * sum * velo;
        }
        case COAG_PROB: {
            double velo = rms_v(u, v);
            double term1 = std::sqrt(3.0 / (2.0 * velo * velo));
            double erfResult = std::erf(term1 * VELO_B);
            double term2 = VELO_B * std::sqrt(6.0 / (M_PI * velo * velo));
            double expResult = std::exp(-3.0 * VELO_B * VELO_B / (2.0 * velo * velo));
            double sum = std::cbrt(u) + std::cbrt(v);
            double sigma_v = M_PI * sum * sum * velo;
            return (erfResult - term2 * expResult) * sigma_v;
        }
        default: throw std::invalid_argument("Unknown coagulation kernel type");
    }
}

/**
 * @brief Fragmentation kernel function
 * @param u First mass
 * @param v Second mass
 * @param fragmentationKernel Selected kernel type
 * @return Kernel value
 */
double fragkernel(double u, double v, FragmentationKernel kernel) {
    switch (kernel) {
        case FRAG_CONSTANT: return 1.0;
        case FRAG_ADDITIVE: return u + v;
        case FRAG_MULTIPLICATIVE: return u * v;
        case FRAG_BALLISTIC: {
            double sum = std::cbrt(u) + std::cbrt(v);
            return M_PI * sum * sum;
        }
        case FRAG_DETER: {
            double velo = rms_v(u, v);
            if (velo >= VELO_F) {
                double sum = std::cbrt(u) + std::cbrt(v);
                return M_PI * sum * sum * velo;
            }
            return 0.0;
        }
        case FRAG_PROB: {
            double velo = rms_v(u, v);
            double term1 = std::sqrt(3.0 / (2.0 * velo * velo));
            double erfResult = std::erf(term1 * VELO_F);
            double term2 = VELO_F * std::sqrt(6.0 / (M_PI * velo * velo));
            double expResult = std::exp(-3.0 * VELO_F * VELO_F / (2.0 * velo * velo));
            double sum = std::cbrt(u) + std::cbrt(v);
            double sigma_v = M_PI * sum * sum * velo;
            return (1.0 - erfResult + term2 * expResult) * sigma_v;
        }
        default: throw std::invalid_argument("Unknown fragmentation kernel type");
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
    
    if (ratio <= 10.0) return total_mass; // Fully destructive
    
    double base_factor = total_mass / (1.0 + ratio);
    
    if (ratio <= 12.0) { // Transition: Destructive -> Cratering
        // f1(q) = 6.5 + 4.5 * cos((q - 10) / 2 * pi)
        double g = 6.5 + 4.5 * std::cos((ratio - 10.0) / 2.0 * M_PI);
        return base_factor * g;
    } else if (ratio <= 15.0) { // Cratering
        return base_factor * 2.0;
    } else if (ratio <= 50.0) { // Transition: Cratering -> Mass Transfer
        // f2(q) = 1.45 + 0.55 * cos((q - 15) / 35 * pi)
        double f = 1.45 + 0.55 * std::cos((ratio - 15.0) / 35.0 * M_PI);
        return base_factor * f;
    } else { // Mass Transfer
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
    if (x < mleft) integral += mleft;
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
    if (x >= mleft) integral += mleft;
    double upper = std::min(mfrag, x);
    integral += mfrag * (std::pow(upper, alpha_p2) - xmin_pow) / (std::pow(mfrag, alpha_p2) - xmin_pow);
    return integral;
}
