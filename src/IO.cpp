#include "IO.hpp"
#include "Globals.hpp"
#include <fstream>
#include <iostream>

// This file contains functions for reading the input file and recording the results to an output file.

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
        std::size_t delimiterPos = line.find('=');
        if (delimiterPos != std::string::npos) {
            std::string key = line.substr(0, delimiterPos);
            std::string value = line.substr(delimiterPos + 1);
            
            // Trim leading and trailing whitespace from key and value
            key.erase(0, key.find_first_not_of(' '));
            key.erase(key.find_last_not_of(' ') + 1);
            value.erase(0, value.find_first_not_of(' '));
            value.erase(value.find_last_not_of(' ') + 1);

            if (key == "NUM_BINS") NUM_BINS = std::stoi(value);
            else if (key == "POLY_ORDER") POLY_ORDER = std::stoi(value);
            else if (key == "MAX_TIME") MAX_TIME = std::stod(value);
            else if (key == "X_MIN") X_MIN = std::stod(value);
            else if (key == "X_MAX") X_MAX = std::stod(value);
            else if (key == "CFL") CFL = std::stod(value);
            else if (key == "VELO_F") VELO_F = std::stod(value);
            else if (key == "VELO_B") VELO_B = std::stod(value);
            else if (key == "ALPHA") ALPHA = std::stod(value);
            else if (key == "QUAD") QUAD = std::stoi(value);
            else if (key == "flux_mode") flux_mode = value;
            else if (key == "COAGULATION_KERNEL") {
                if (value == "constant") coagulationKernel = COAG_CONSTANT;
                else if (value == "additive") coagulationKernel = COAG_ADDITIVE;
                else if (value == "multiplicative") coagulationKernel = COAG_MULTIPLICATIVE;
                else if (value == "ballistic") coagulationKernel = COAG_BALLISTIC;
                else if (value == "determined") coagulationKernel = COAG_DETER;
                else if (value == "probabilistic") coagulationKernel = COAG_PROB;
                else {std::cerr << "Error: Unknown COAGULATION_KERNEL value: " << value << std::endl;
                      exit(1);}
                // Add parsing for others
            }
            else if (key == "FRAGMENTATION_KERNEL") {
                if (value == "constant") fragmentationKernel = FRAG_CONSTANT;
                else if (value == "additive") fragmentationKernel = FRAG_ADDITIVE;
                else if (value == "multiplicative") fragmentationKernel = FRAG_MULTIPLICATIVE;
                else if (value == "ballistic") fragmentationKernel = FRAG_BALLISTIC;
                else if (value == "determined") fragmentationKernel = FRAG_DETER;
                else if (value == "probabilistic") fragmentationKernel = FRAG_PROB;
                else {std::cerr << "Error: Unknown FRAGMENTATION_KERNEL value: " << value << std::endl;
                      exit(1);}
                // Add parsing for others
            }
        }
    }
    inputFile.close();
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
void record(const std::string& filename, std::vector<std::vector<double>>& g, std::vector<double>& cell, std::vector<double>& h, std::vector<double>& x, int step, double T) {
    std::ofstream fp;
    // If step == 0, overwrite the file; else append to it
    if (step == 0) fp.open(filename, std::ios::trunc); 
    else fp.open(filename, std::ios::app);

    fp << "Step: " << step << "\nTime: " << T << "\n";
    // If it's the first step (step == 0), save cell, h, and x
    if (step == 0) {
        for (int i = 0; i < NUM_BINS + 1; i++) fp << cell[i] << "\n";
        for (int i = 0; i < NUM_BINS; i++) fp << h[i] << "\n";
        for (int i = 0; i < NUM_BINS; i++) fp << x[i] << "\n";
    }
    for (int i = 0; i < NUM_BINS; i++) {
        for (int j = 0; j < POLY_ORDER + 1; j++) fp << g[i][j] << " ";
        fp << "\n";
    }
    fp.close();
}
