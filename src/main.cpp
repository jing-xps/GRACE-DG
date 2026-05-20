/// @file main.cpp
/// @brief Main entry point for the GRACE-DG numerical solver.

#include <iostream>
#include <vector>
#include <chrono>
#include "Globals.hpp"
#include "IO.hpp"
#include "MathUtils.hpp"
#include "Physics.hpp"
#include "DGSolver.hpp"

/**
 * @brief Main simulation driver
 * @return Exit status
 */
int main() {
    TimeProfile timer;
    readInputFile("input.txt");
   
    // Print configuration
    std::cout << "--- GRACE-DG Setup ---\n"
              << "NUM_BINS: " << NUM_BINS << "\n"
              << "POLY_ORDER: " << POLY_ORDER << "\n"
              << "QUAD: " << QUAD << "\n"
              << "MAX_TIME: " << MAX_TIME << "\n"
              << "X_MIN: " << X_MIN << "\n"
              << "X_MAX: " << X_MAX << "\n"
              << "CFL: " << CFL << "\n"
              << "VELO_F: " << VELO_F << "\n"
              << "VELO_B: " << VELO_B << "\n"
              << "ALPHA: " << ALPHA << "\n"
              << "Flux Mode: " << flux_mode << std::endl;
    
    // Initialize quadrature points, weights, and precomputed constants
    init_quadrature_globals();
    init_precomputed_constants(); 
   
    // Initialize laqp and Dlaqp
    std::vector<std::vector<double>> laqp = init_laqp();
    std::vector<std::vector<double>> Dlaqp = init_Dlaqp();
    
    // Initialize mass bins, cell boundaries, and related variables
    std::vector<double> h(NUM_BINS), x(NUM_BINS), cell(NUM_BINS + 1);
    initialize_mass_bin(h, x, cell);
    
    // Initialize density and related variables
    std::vector<std::vector<double>> g(NUM_BINS, std::vector<double>(POLY_ORDER + 1));
    std::vector<double> g_bar(NUM_BINS);
    std::vector<double> gamma(NUM_BINS);
    initialize_mass_density(g, h, x, cell);
    cal_gamma(g, h, x, cell, g_bar, gamma);
    std::vector<std::vector<double>> p(NUM_BINS, std::vector<double>(POLY_ORDER + 1));
    cal_p(g, gamma, p);

    
    int step = 0;
    double T = 0;
    std::vector<double> flux(NUM_BINS + 1);
    // Record initial state at step 0 (overwrite the file)
    record("output.txt", p, cell, h, x, step, T);
   
    // Main loop for simulation
    while (T < MAX_TIME) {
        auto loop_start = std::chrono::high_resolution_clock::now();
        
        // Calculate time step
        auto t1 = std::chrono::high_resolution_clock::now();  
        double dt = cal_time_step(p, h, x, cell, flux, laqp);
        timer.total_cal_time_step += std::chrono::duration<double>(std::chrono::high_resolution_clock::now() - t1).count();
        
        if (T + dt > MAX_TIME) dt = MAX_TIME - T;
        
        bool success = false;
        while (!success) {
            // Evolve system
            auto t3 = std::chrono::high_resolution_clock::now();
            success = evolve(g, p, h, x, cell, flux, dt, laqp, Dlaqp);
            timer.total_evolve += std::chrono::duration<double>(std::chrono::high_resolution_clock::now() - t3).count();
        }
       
        // Recalculate gamma and update p
        auto t5 = std::chrono::high_resolution_clock::now();
        cal_gamma(g, h, x, cell, g_bar, gamma);
        timer.total_cal_gamma += std::chrono::duration<double>(std::chrono::high_resolution_clock::now() - t5).count();
        

        auto t7 = std::chrono::high_resolution_clock::now();
        cal_p(g, gamma, p);
        timer.total_cal_p += std::chrono::duration<double>(std::chrono::high_resolution_clock::now() - t7).count();
        
        // Update time and step count
        T += dt;
        step++;
        
        // Record current state (append to the file)
        auto t9 = std::chrono::high_resolution_clock::now();
        record("output.txt", p, cell, h, x, step, T);
        timer.total_record += std::chrono::duration<double>(std::chrono::high_resolution_clock::now() - t9).count();

        // Update total loop time
        timer.total_loop += std::chrono::duration<double>(std::chrono::high_resolution_clock::now() - loop_start).count();
        
        std::cout << "Step: " << step << " | Time: " << T << " | dt: " << dt << std::endl;
    }
    
    // Print performance profile
    std::cout << "\n--- Performance Profile (Step " << step << ") ---\n"
              << "cal_time_step: " << timer.total_cal_time_step << " s\n"
              << "evolve: " << timer.total_evolve << " s\n"
              << "cal_gamma: " << timer.total_cal_gamma << " s\n"
              << "cal_p: " << timer.total_cal_p << " s\n"
              << "record: " << timer.total_record << " s\n"
              << "Total loop: " << timer.total_loop << " s\n";
              
    return 0;
}
