#pragma once
#include <vector>

void initialize_mass_bin(std::vector<double>& h, std::vector<double>& x, std::vector<double>& cell);
void initialize_mass_density(std::vector<std::vector<double>>& g, std::vector<double>& h, std::vector<double>& x, std::vector<double>& cell);
void cal_cell_average(std::vector<std::vector<double>>& g, std::vector<double>& g_bar);
void cal_gamma(std::vector<std::vector<double>>& g, std::vector<double>& h, std::vector<double>& x, std::vector<double>& cell, std::vector<double>& g_bar, std::vector<double>& gamma);
void cal_p(std::vector<std::vector<double>>& g, std::vector<double>& gamma, std::vector<std::vector<double>>& p);

// Flux calculations (non-conservative and conservative)
double cal_flux_coag_nc(std::vector<std::vector<double>>& p, std::vector<double>& h, std::vector<double>& x, std::vector<double>& cell, int j, const std::vector<std::vector<double>>& laqp);
double cal_flux_integral_coag_nc(std::vector<std::vector<double>>& p, std::vector<double>& h, std::vector<double>& x, std::vector<double>& cell, int j, int i, const std::vector<std::vector<double>>& Dlaqp);
double cal_flux_frag_nc(std::vector<std::vector<double>>& p, std::vector<double>& h, std::vector<double>& x, std::vector<double>& cell, int j, const std::vector<std::vector<double>>& laqp);
double cal_flux_integral_frag_nc(std::vector<std::vector<double>>& p, std::vector<double>& h, std::vector<double>& x, std::vector<double>& cell, int j, int i, const std::vector<std::vector<double>>& laqp);

double cal_flux_coag_c(std::vector<std::vector<double>>& p, std::vector<double>& h, std::vector<double>& x, std::vector<double>& cell, int j, const std::vector<std::vector<double>>& laqp);
double cal_flux_integral_coag_c(std::vector<std::vector<double>>& p, std::vector<double>& h, std::vector<double>& x, std::vector<double>& cell, int j, int i, const std::vector<std::vector<double>>& Dlaqp);
double cal_flux_frag_c(std::vector<std::vector<double>>& p, std::vector<double>& h, std::vector<double>& x, std::vector<double>& cell, int j, const std::vector<std::vector<double>>& laqp);
double cal_flux_integral_frag_c(std::vector<std::vector<double>>& p, std::vector<double>& h, std::vector<double>& x, std::vector<double>& cell, int j, int i, const std::vector<std::vector<double>>& laqp);

double cal_time_step(std::vector<std::vector<double>>& p, std::vector<double>& h, std::vector<double>& x, std::vector<double>& cell, std::vector<double>& flux, const std::vector<std::vector<double>>& laqp);
bool evolve(std::vector<std::vector<double>>& g, std::vector<std::vector<double>>& p, std::vector<double>& h, std::vector<double>& x, std::vector<double>& cell, std::vector<double>& flux, double& dt, const std::vector<std::vector<double>>& laqp, const std::vector<std::vector<double>>& Dlaqp);
