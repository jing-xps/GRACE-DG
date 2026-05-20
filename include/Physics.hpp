#pragma once
#include "Globals.hpp"

/// @brief Initial mass density function
double g_0(double x);

/// @brief Root mean square relative velocity function
double rms_v(double u, double v);

/// @brief Coagulation kernel function
double coagkernel(double u, double v, CoagulationKernel kernel);

/// @brief Fragmentation kernel function
double fragkernel(double u, double v, FragmentationKernel kernel);

/// @brief Calculate total fragment mass after collision
double m_frag(double u, double v);

/// @brief Helper function for fragmentation integral (type 1)
double inte1(double u, double v, double x);

/// @brief Helper function for fragmentation integral (type 2)
double inte2(double u, double v, double x);
