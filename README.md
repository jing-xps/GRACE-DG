# GRACE-DG (GRAin collisional Evolution - Discontinuous Galerkin)

GRACE-DG is a high-order numerical solver for the Smoluchowski coagulation-fragmentation equation. Utilizing the Discontinuous Galerkin (DG) method, it is designed for accurate and efficient modeling of dust grain collisional evolution, making it highly applicable to environments such as protoplanetary disks.

## Features

- **High-Order Accuracy:** Implements the Discontinuous Galerkin method up to the 4th polynomial order.
- **Flexible Physics:** Supports various coagulation and fragmentation kernels (Constant, Additive, Multiplicative, Ballistic, Deterministic, Probabilistic).
- **Adaptive Time-stepping:** Automatically adjusts the time step based on flux conditions (CFL constraint).
- **Modular Design:** Easy to extend with custom fragmentation distributions or collision velocities.

## Build Instructions

This project uses CMake. To build the solver, run the following commands in the project root directory:

```bash
mkdir build
cd build
cmake ..
make

cd ../examples/coag0
cp ../../build/GRACE_DG .
./GRACE_DG
```
