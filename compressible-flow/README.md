# 1. Godunov Scheme
## 1.1. Definition
- Piecewise Constant Approximation: Assuming that the solution at each time step is piecewise constant, dividing the domain into cells.
- Approximate Riemann Solver: At each cell interface, the scheme solves the Riemann problem.
- Finite Volume Formulation: The Godunov scheme is typically implemented in a finite volume framework, where the fluxes across the boundaries of each control volume are computed to update the cell-averaged values.
- First-Order Accuracy
- Upwind Nature
  
