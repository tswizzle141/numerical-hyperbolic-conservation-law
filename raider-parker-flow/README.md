# 1. Raider–Parker Problem
This problem studies the behavior of one-dimensional, inviscid compressible flows. Its setup often resembles Sod's Shock Tube Problem, where a sudden discontinuity in density and pressure propagates through the domain, generating waves (shocks, rarefactions, or contact discontinuities). Mathematically, using the Euler equations for compressible flow:
- Mass Conservation (Continuity Equation):
$$\frac{\partial \rho}{\partial t} + \frac{\partial (\rho u)}{\partial x}=0$$
- Momentum Conservation:
$$\frac{\partial (\rho u)}{\partial t} + \frac{\partial (\rho u^2+p)}{\partial x}=0$$
- Energy Conservation:
$$\frac{\partial E}{\partial t} +\frac{\partial [(E+p)u]}{\partial x}=0$$

with $\rho$ is density, $u$ is velocity, $p$ is pressure, $E=\rho e+\frac{1}{2}\rho u^2$ is total energy per unit volumn, $e=\frac{p}{(\gamma-1)\rho}$ is internal energy per unit mass, $\gamma$ is ratio of specific heat.

# 2. Lax-Wendroff Scheme
Lax-Wendroff scheme is a second-order numerical method solving hyperbolic PDEs, combining Taylor expansion with FDM. For conserved variable $U(x,t)$, this scheme uses a second-order Taylor expansion in time:

$$U_i^{n+1}=U_i^n + \Delta t \frac{\partial U}{\partial t} + \frac{\Delta t^2}{2}\frac{\partial^2 U}{\partial t^2}$$

# 3. Solve Raider-Parker problem step-by-step
* Step 1: Define initial conditions for density $\rho$, velocity $u$, pressure $p$ in the domain. In this problem domain is often divided into 2 regions:
  - Left region: High density and pressure
  - Right region: Low density and pressure
* Step 2: Discretization
Divide spatial domain with step size $\Delta x$, choose time step $\Delta t$ based on CFL condition: $\Delta t=\text{CFL} \cdot \frac{\Delta x}{\max(|u|+c)}$ where $c=\sqrt{\frac{\gamma p}{\rho}}$.
* Step 3: Compute flux
Use Euler function to define flux vector $F(U)$:
$$F(U) = \begin{bmatrix} \rho u\\
        \rho u^2+p\\
        (E+p)u
        \end{bmatrix}$$

* Step 4: 
