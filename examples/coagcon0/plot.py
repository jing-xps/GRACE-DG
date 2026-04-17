import numpy as np
import matplotlib.pyplot as plt
import os
from math import *
from scipy.special import iv  # I1 is the modified Bessel function of the first kind

def legendre(i, x):
    if i == 0:
        return 1
    if i == 1:
        return x
    Pim1 = 1
    Pi = x
    for n in range(1, i):
        Pip1 = ((2 * n + 1) * x * Pi - n * Pim1) / (n + 1)
        Pim1 = Pi
        Pi = Pip1
    return Pi


# Define the analytical solution; for constant kernel and g0 = x * exp(-x);
def analytic_g(x, t, kernel):
    if kernel == 'constant':
        return (4 * x) / (2 + t)**2 * np.exp(-(1 - t / (2 + t)) * x)
    elif kernel == 'additive':
        T = 1 - np.exp(-t)  # Define T
        term1 = (1 - T) * np.exp(-x * (1 + T))  # First part of the numerator
        term2 = T**(1/2)  # Denominator
        bessel_term = iv(1, 2 * x * T**(1/2))  # Modified Bessel function of the first kind
        return term1 / term2 * bessel_term
    elif kernel == 'multiplicative':
        # Piecewise definition of T
        if t <= 1:
            T = 1 + t
        else:
            T = 2 * t**(1/2)
        term1 = np.exp(-T * x)  # Exponential term
        bessel_term = iv(1, 2 * x * t**(1/2))  # Modified Bessel function of the first kind
        term2 = x**2 * t**(1/2)  # Denominator
        return x* term1 * bessel_term / term2
    else:
        print('Unknown coagulation kernel type')

# Function to read parameters from input.txt
def read_parameters(filename):
    params = {}
    with open(filename, 'r') as file:
        for line in file:
            line = line.strip()
            if line:  # Skip empty lines
                if "COAGULATION_KERNEL" not in line and "FRAGMENTATION_KERNE" not in line and "flux_mode" not in line:
                    key, value = line.split('=')
                    params[key.strip()] = float(value.strip())
                else:
                    key, value = line.split('=')
                    params[key.strip()] = str(value.strip())
    return params

# Read parameters from input.txt
params = read_parameters('input.txt')
print(params)

NUM_BINS = int(params['NUM_BINS'])
POLY_ORDER = int(params['POLY_ORDER'])
MAX_TIME = params['MAX_TIME']
X_MIN = params['X_MIN']
X_MAX = params['X_MAX']
kernel = params['COAGULATION_KERNEL']

# Ensure output directory exists
if not os.path.exists('./pic'):
    os.makedirs('./pic')

# Open the output file and read all lines
with open('output.txt', 'r') as file:
    lines = file.readlines()

# Initialize variables to hold cell, h, x values
cell, h, x = None, None, None

# Loop through each step (starting from step 0)
step = 0
T = 0

while True:
    if not lines:
        break
    # Read the step number
    step_line = lines.pop(0)
    step = int(step_line.split(": ")[1])

    # Read the time T for the current step
    time_line = lines.pop(0)
    T = float(time_line.split(": ")[1])  # Extract the time value

    # Read cell, h, x, and g values for step 0
    if step == 0:
        cell = np.array([float(lines.pop(0)) for _ in range(NUM_BINS + 1)])
        h = np.array([float(lines.pop(0)) for _ in range(NUM_BINS)])
        x = np.array([float(lines.pop(0)) for _ in range(NUM_BINS)])

    # Read g values for the current step
    print("step: ", step)
    print("time ", T)
    g = np.zeros((NUM_BINS, POLY_ORDER + 1))
    for i in range(NUM_BINS):
        g[i] = np.array([float(val) for val in lines.pop(0).split()])

    # Plot each bin's function
    plt.figure(figsize=(8, 6))
    for c in cell:
        plt.axvline(x=c, color='gray', linestyle='--', linewidth=1)
    
    for i in range(NUM_BINS):
        x_values = np.linspace(cell[i], cell[i+1], 100)
        y_values = np.zeros_like(x_values)
        
        # Sum up the Legendre polynomials for the current cell
        for k in range(POLY_ORDER + 1):
            y_values += g[i, k] * legendre(k, 2/h[i] * (x_values - x[i]))
        if i == 0:
            plt.plot(x_values, y_values, color='r', linewidth=2, label=f'Numerical solution (k={POLY_ORDER})')
            # Plot for the current cell
        else:
            plt.plot(x_values, y_values, color='r', linewidth=2)
    
    # Plot analytical g(x) for the current step
    x_func = np.logspace(log10(X_MIN), log10(X_MAX), 500)
    y_func = analytic_g(x_func, T, kernel)  # Use the extracted time T
    plt.plot(x_func, y_func, color='black', linestyle='--', label='Analytic solution', linewidth=2)

    # Set the logarithmic scale for the x-axis and labels
    plt.xscale('log')
    plt.xlabel(r'mass $x$', fontsize=15)
    plt.ylabel(r'mass density $g(x,\tau)$', fontsize=15)

    # Add a legend and save the figure
    # plt.legend([f'k={POLY_ORDER}'], loc='best')
    plt.title(f"Step {step}, Time {T:.0f}")
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    plt.legend(loc='best', fontsize=15)
    plt.xlim(X_MIN,X_MAX)
    # Save the figure as a PNG file
    plt.savefig(f'./pic/{step}.png', dpi=300)

    # Close the plot to free memory
    plt.close()



print("Plots have been saved in the ./pic/ directory.")

