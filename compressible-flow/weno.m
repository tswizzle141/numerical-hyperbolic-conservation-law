import numpy as np
import matplotlib.pyplot as plt

def compute_flux(U, gamma):
    """Compute flux vector from conserved variables."""
    rho = U[0]
    u = U[1] / rho
    E = U[2]
    p = (gamma - 1) * (E - 0.5 * rho * u**2)
    return np.array([rho * u, rho * u**2 + p, (E + p) * u])

def weno_reconstruct(U, direction, epsilon=1e-6):
    """Perform WENO5 reconstruction."""
    def weights(beta):
        alpha = [1 / (epsilon + b)**2 for b in beta]
        return [a / sum(alpha) for a in alpha]

    def smoothness_indicators(stencil):
        beta = [
            (13 / 12) * (stencil[0] - 2 * stencil[1] + stencil[2])**2 +
            (1 / 4) * (stencil[0] - 4 * stencil[1] + 3 * stencil[2])**2,
            (13 / 12) * (stencil[1] - 2 * stencil[2] + stencil[3])**2 +
            (1 / 4) * (stencil[1] - stencil[3])**2,
            (13 / 12) * (stencil[2] - 2 * stencil[3] + stencil[4])**2 +
            (1 / 4) * (3 * stencil[2] - 4 * stencil[3] + stencil[4])**2
        ]
        return beta

    n = len(U)
    UL = np.zeros(n - 4)
    UR = np.zeros(n - 4)

    for i in range(2, n - 2):
        stencil = U[i - 2:i + 3]
        beta = smoothness_indicators(stencil)
        w = weights(beta)

        if direction == "left":
            UL[i - 2] = w[0] * (2 * stencil[0] - stencil[1] + stencil[2]) / 6 + \
                        w[1] * (stencil[1] + stencil[2] + stencil[3]) / 6 + \
                        w[2] * (-stencil[2] + 5 * stencil[3] - stencil[4]) / 6
        elif direction == "right":
            UR[i - 2] = w[0] * (2 * stencil[4] - stencil[3] + stencil[2]) / 6 + \
                        w[1] * (stencil[3] + stencil[2] + stencil[1]) / 6 + \
                        w[2] * (-stencil[2] + 5 * stencil[1] - stencil[0]) / 6

    return UL if direction == "left" else UR

def hllc_flux(UL, UR, gamma):
    """Compute flux at the interface using the HLLC solver."""
    rhoL, uL, pL = UL
    rhoR, uR, pR = UR

    # Compute wave speeds
    cL = np.sqrt(gamma * pL / rhoL)
    cR = np.sqrt(gamma * pR / rhoR)
    SL = min(uL - cL, uR - cR)
    SR = max(uL + cL, uR + cR)

    if SL >= 0:
        return compute_flux(UL, gamma)
    elif SR <= 0:
        return compute_flux(UR, gamma)
    else:
        # Compute flux in the star region
        S_star = (pR - pL + rhoL * uL * (SL - uL) - rhoR * uR * (SR - uR)) / \
                 (rhoL * (SL - uL) - rhoR * (SR - uR))
        FL = compute_flux(UL, gamma)
        FR = compute_flux(UR, gamma)
        return (SR * FL - SL * FR + SL * SR * (UR - UL)) / (SR - SL)

def sod_shock_tube(gamma, nx, dx, dt, t_max, CFL):
    """Solve Sod's Shock Tube problem using WENO scheme."""
    # Initial conditions
    U = np.zeros((3, nx))
    rho = np.where(np.arange(nx) < nx // 2, 1.0, 0.125)
    p = np.where(np.arange(nx) < nx // 2, 1.0, 0.1)
    u = np.zeros(nx)

    E = p / (gamma - 1) + 0.5 * rho * u**2
    U[0, :] = rho
    U[1, :] = rho * u
    U[2, :] = E

    t = 0
    while t < t_max:
        # Compute time step
        a = np.sqrt(gamma * p / rho)  # Speed of sound
        dt = CFL * dx / np.max(np.abs(u) + a)
        if t + dt > t_max:
            dt = t_max - t

        # Add ghost cells for WENO reconstruction
        U_ext = np.pad(U, ((0, 0), (2, 2)), mode="edge")

        # Compute fluxes
        F = np.zeros_like(U)
        for i in range(2, nx + 2):
            UL = weno_reconstruct(U_ext[:, i - 2:i + 3], "left")
            UR = weno_reconstruct(U_ext[:, i - 2:i + 3], "right")
            F[:, i - 2] = hllc_flux(UL, UR, gamma)

        # Update conserved variables
        U[:, 1:-1] -= dt / dx * (F[:, 1:] - F[:, :-1])

        # Update primitive variables
        rho = U[0, :]
        u = U[1, :] / rho
        E = U[2, :]
        p = (gamma - 1) * (E - 0.5 * rho * u**2)

        t += dt

    return rho, u, p

# Parameters
gamma = 1.4
nx = 200
dx = 1.0 / nx
t_max = 0.2
CFL = 0.8

# Solve the problem
rho, u, p = sod_shock_tube(gamma, nx, dx, 0.001, t_max, CFL)

# Plot the results
x = np.linspace(0, 1, nx)
plt.figure(figsize=(12, 8))
plt.subplot(311)
plt.plot(x, rho, label="Density")
plt.legend()
plt.subplot(312)
plt.plot(x, u, label="Velocity")
plt.legend()
plt.subplot(313)
plt.plot(x, p, label="Pressure")
plt.legend()
plt.xlabel("x")
plt.tight_layout()
plt.show()
