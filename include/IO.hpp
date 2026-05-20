#pragma once
#include <string>
#include <vector>

void readInputFile(const std::string& filename);
void record(const std::string& filename, std::vector<std::vector<double>>& g, std::vector<double>& cell, std::vector<double>& h, std::vector<double>& x, int step, double T);
