# 1. Godunov Scheme
## 1.1. Definition
- Piecewise Constant Approximation: Assuming that the solution at each time step is piecewise constant, dividing the domain into cells.
- Approximate Riemann Solver: At each cell interface, the scheme solves the Riemann problem.
- Finite Volume Formulation: The Godunov scheme is typically implemented in a finite volume framework, where the fluxes across the boundaries of each control volume are computed to update the cell-averaged values.
- First-Order Accuracy
- Upwind Nature
## 1.2. Properties
- Governed by conservation law
- Capture shock waves and discontinuities without introducing spurious oscillations
- Stable under CFL condition
- 1st-order Godunov scheme is monotone

# 2. Problem on Compressible Flows
## 2.1. Problem Statement (Sod's Shock Tube Problem - 1D Compressible Flow)
$$\frac{\partial}{\partial t} \begin{bmatrix} \rho\\
\rho u\\
E\\
\end{bmatrix} + \frac{\partial}{\partial x} \begin{bmatrix} \rho u\\
\rho u^2+p\\
u(E+p)\\
\end{bmatrix}=0$$

where $\rho$ is density, $u$ is velocity, $p$ is pressure, $E=\frac{p}{\gamma-1}+\frac{1}{2}\rho u^2$ is total energy per unit volume, with $\gamma$ is specific ratio.
## 2.2. Step-by-step solution
* Step 1: Domain discretization
Divide domain into $N$ cells with spatial step size $\Delta x$. Let cell averages of conserved quantities $(U)$ at time $t$ be:
$$U_i^n = \frac{1}{\Delta x} \displaystyle\int_{x_{i-\frac{1}{2}}}^{x_{i+\frac{1}{2}}} U(x,t^n)dx$$ where $U=[\rho, \rho u, E]^T$ is conserved variables. 
* Step 2: Solve Riemann problem
At each cell interface $x_{i+\frac{1}{2}}$, solve Riemann problem to determine the flux $F_{i+\frac{1}{2}}^n$:

$$F=\begin{cases} \rho u \\
                  \rho u^2+p \\
                  u(E+p) \end{cases}$$

The Riemann solver uses $U_L, U_R$ states at the interface.
* Step 3: Update the Conserved Quantities $$U_i^{n+1}=U_i^n-\frac{\Delta t}{\Delta x}\Big(F_{i+\frac{1}{2}}^n-F_{i-\frac{1}{2}}^n \Big)$$
* 

$$\rho(x,0) = \begin{cases} 1 \ \ \ x < 0.5 \\
                          0.125 \ \ \ x \geq 0.5
\end{cases}$$
$$u(x,0) = 0$$
$$p(x,0) = \begin{cases} 1 \ \ \ x < 0.5 \\
                        0.1 \ \ \ x \geq 0.5
\end{cases}$$

# 3. MUSCL Scheme
# 4. WENO Scheme
