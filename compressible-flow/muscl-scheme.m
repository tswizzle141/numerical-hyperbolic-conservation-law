%parameters for Sod's Shock Tube
Nx = 200; %number of grid points
x = linspace(0, 1, Nx); %spatial grid
dx = x(2) - x(1); %step size
gamma = 1.4; %heat ratio for air
CFL = 0.9; %CFL condition
t_final = 0.2;

%initial conditions (Left and Right states)
rho = ones(1, Nx); %density
u = zeros(1, Nx); %velocity
p = ones(1, Nx); %pressure

%left region (x<0.5)
rho(x < 0.5) = 1.0;
u(x < 0.5) = 0.0;
p(x < 0.5) = 1.0;

%right region (x>=0.5)
rho(x >= 0.5) = 0.125;
u(x >= 0.5) = 0.0;
p(x >= 0.5) = 0.1;

%conservative variables
E = p / (gamma - 1) + 0.5 * rho .* u.^2; %total energy
U = [rho; rho .* u; E];                  %conservative variables

%time-stepping loop
t = 0;
while t < t_final
    %calculate time step based on CFL condition
    a = sqrt(gamma * p ./ rho); %sound speed
    dt = CFL * dx / max(abs(u) + a);
    if t + dt > t_final
        dt = t_final - t; %adjust for final time step
    end

    %MUSCL reconstruction: Compute left and right states at interfaces
    UL = zeros(size(U));
    UR = zeros(size(U));
    for i = 2:Nx-1
        for var = 1:3 %loop over conservative variables
            %compute slopes using minmod limiter
            slope_left = U(var, i) - U(var, i-1);
            slope_right = U(var, i+1) - U(var, i);
            slope = minmod(slope_left, slope_right);

            %reconstruct left and right states
            UL(var, i) = U(var, i) - 0.5 * slope; %left state
            UR(var, i) = U(var, i) + 0.5 * slope; %right state
        end
    end

    %apply boundary conditions (ghost cells)
    UL(:, 1) = U(:, 1);
    UR(:, Nx) = U(:, Nx);

    %compute fluxes using HLLC Riemann solver
    F = zeros(size(U));
    for i = 2:Nx
        F(:, i-1) = HLLC_flux(UL(:, i-1), UR(:, i), gamma);
    end

    %update conservative variables
    for i = 2:Nx-1
        U(:, i) = U(:, i) - dt / dx * (F(:, i) - F(:, i-1));
    end

    %update primitive variables
    rho = U(1, :);
    u = U(2, :) ./ rho;
    E = U(3, :);
    p = (gamma - 1) * (E - 0.5 * rho .* u.^2);

    t = t + dt;
end

figure;
subplot(3, 1, 1);
plot(x, rho, 'LineWidth', 2);
title('Density');
xlabel('x');
ylabel('\rho');
grid on;

subplot(3, 1, 2);
plot(x, u, 'LineWidth', 2);
title('Velocity');
xlabel('x');
ylabel('u');
grid on;

subplot(3, 1, 3);
plot(x, p, 'LineWidth', 2);
title('Pressure');
xlabel('x');
ylabel('p');
grid on;

%minmod limiter function
function result = minmod(a, b)
    result = 0.5 * (sign(a) + sign(b)) .* min(abs(a), abs(b));
end

%HLLC flux function
function F = HLLC_flux(UL, UR, gamma)
    %extract left and right states
    rhoL = UL(1); uL = UL(2) / rhoL; EL = UL(3);
    rhoR = UR(1); uR = UR(2) / rhoR; ER = UR(3);

    %compute pressure and speed of sound
    pL = (gamma - 1) * (EL - 0.5 * rhoL * uL^2);
    pR = (gamma - 1) * (ER - 0.5 * rhoR * uR^2);
    aL = sqrt(gamma * pL / rhoL);
    aR = sqrt(gamma * pR / rhoR);

    %compute wave speeds
    SL = min(uL - aL, uR - aR);
    SR = max(uL + aL, uR + aR);
    S_star = (pR - pL + rhoL * uL * (SL - uL) - rhoR * uR * (SR - uR)) / ...
             (rhoL * (SL - uL) - rhoR * (SR - uR));

    %compute fluxes for left, right, and star regions
    FL = [rhoL * uL; rhoL * uL^2 + pL; uL * (EL + pL)];
    FR = [rhoR * uR; rhoR * uR^2 + pR; uR * (ER + pR)];
    if 0 <= SL
        F = FL;
    elseif SL <= 0 && 0 <= S_star
        U_star_L = [rhoL * (SL - uL) / (SL - S_star); ...
                    rhoL * (SL - uL) / (SL - S_star) * S_star; ...
                    (SL - uL) / (SL - S_star) * (EL + (S_star - uL) * (S_star + pL / (rhoL * (SL - uL))))];
        F = FL + SL * (U_star_L - UL);
    elseif S_star <= 0 && 0 <= SR
        U_star_R = [rhoR * (SR - uR) / (SR - S_star); ...
                    rhoR * (SR - uR) / (SR - S_star) * S_star; ...
                    (SR - uR) / (SR - S_star) * (ER + (S_star - uR) * (S_star + pR / (rhoR * (SR - uR))))];
        F = FR + SR * (U_star_R - UR);
    else
        F = FR;
    end
end
