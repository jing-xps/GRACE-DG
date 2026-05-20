#pragma once
#include <string>
#include <vector>

/// @brief Structure for tracking computation time of different components
struct TimeProfile {
    double total_cal_time_step = 0;
    double total_evolve = 0;
    double total_cal_gamma = 0;
    double total_cal_p = 0;
    double total_record = 0;
    double total_loop = 0;
};

// Coagulation kernel types
enum CoagulationKernel {
    COAG_CONSTANT,        ///< Constant kernel: K = 1 
    COAG_ADDITIVE,        ///< Additive kernel: K = u + v
    COAG_MULTIPLICATIVE,  ///< Multiplicative kernel: K = u * v
    COAG_BALLISTIC,       ///< Ballistic kernel: K = pi * (u^(1/3) + v^(1/3))^2
    COAG_DETER,           ///< Deterministic kernel with velocity threshold
    COAG_PROB,            ///< Probabilistic kernel with velocity threshold
};

// Fragmentation kernel types; same as coagulation but for fragmentation processes
enum FragmentationKernel {
    FRAG_CONSTANT,       
    FRAG_ADDITIVE,       
    FRAG_MULTIPLICATIVE, 
    FRAG_BALLISTIC,      
    FRAG_DETER,          
    FRAG_PROB,           
};

// Global Configuration Parameters
extern int NUM_BINS;
extern int POLY_ORDER;
extern int QUAD;
extern double MAX_TIME;
extern double X_MIN;
extern double X_MAX;
extern double CFL;
extern double VELO_F;
extern double VELO_B;
extern double ALPHA;

extern CoagulationKernel coagulationKernel;
extern FragmentationKernel fragmentationKernel;
extern std::string flux_mode;

// Quadrature global variables and precomputed constants
extern const double* Q_POINTS;
extern const double* Q_WEIGHTS;
extern int Q_SIZE;
extern double alpha_p2;
extern double xmin_pow;

/// Initialize quadrature information and constants
void init_quadrature_globals();
void init_precomputed_constants();
