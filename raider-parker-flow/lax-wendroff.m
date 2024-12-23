%Raider-Parker Flow using Lax-Wendroff Scheme
clear; clc; close all;

% Parameters
nx = 200; %number of grid points
x = linspace(0, 1, nx); %spatial domain
dx = x(2) - x(1); %step size
gamma = 1.4; %ratio of specific heats for air
t_max = 0.15;
CFL = 0.8; %CFL number

%initial conditions (Raider-Parker-like initial state)
rho = ones(1, nx); %density
u = zeros(1, nx); %velocity
p = ones(1, nx); %pressure
rho(x > 0.5) = 0.125; %lower density region
p(x > 0.5) = 0.1; %lower pressure region

%derived quantities
E = p / (gamma - 1) + 0.5 * rho .* u.^2; %total energy
U = [rho; rho .* u; E]; %conserved variables

%time stepping
t = 0;
dt = CFL * dx / max(sqrt(gamma * p ./ rho) + abs(u)); %initial time step

while t < t_max
    %adjust time step for final iteration
    if t + dt > t_max
        dt = t_max - t;
    end

    %compute fluxes
    F = flux(U, gamma);

    %Lax-Wendroff Scheme
    U_half = 0.5 * (U(:, 1:end-1) + U(:, 2:end)) - ...
             0.5 * dt / dx * (F(:, 2:end) - F(:, 1:end-1));
    F_half = flux(U_half, gamma);

    U(:, 2:end-1) = U(:, 2:end-1) - dt / dx * (F_half(:, 2:end) - F_half(:, 1:end-1));

    %update primitive variables
    rho = U(1, :);
    u = U(2, :) ./ rho;
    E = U(3, :);
    p = (gamma - 1) * (E - 0.5 * rho .* u.^2);

    %update time
    t = t + dt;
    dt = CFL * dx / max(sqrt(gamma * p ./ rho) + abs(u)); % Recalculate time step
end

figure;
subplot(3, 1, 1);
plot(x, rho, 'b-', 'LineWidth', 1.5);
xlabel('x'); ylabel('Density');
title('Density Distribution'); grid on;

subplot(3, 1, 2);
plot(x, u, 'r-', 'LineWidth', 1.5);
xlabel('x'); ylabel('Velocity');
title('Velocity Distribution'); grid on;

subplot(3, 1, 3);
plot(x, p, 'k-', 'LineWidth', 1.5);
xlabel('x'); ylabel('Pressure');
title('Pressure Distribution'); grid on;

%compute flux
function F = flux(U, gamma)
    rho = U(1, :);
    u = U(2, :) ./ rho;
    E = U(3, :);
    p = (gamma - 1) * (E - 0.5 * rho .* u.^2);
    F = [rho .* u; rho .* u.^2 + p; (E + p) .* u];
end
