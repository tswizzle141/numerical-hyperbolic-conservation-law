%params for Sod's shock tube
Nx = 200; %number of grid points
x = linspace(0, 1, Nx); %spatial grid
dx = x(2) - x(1); %step size
gamma = 1.4; %heat ratio for air
CFL = 0.9; %CFL condition
t_final = 0.2;           

%initial conditions (left and right states)
rho = ones(1, Nx); %density
u = zeros(1, Nx); %velocity
p = ones(1, Nx); %pressure

%left region (x < 0.5)
rho(x < 0.5) = 1.0;
u(x < 0.5) = 0.0;
p(x < 0.5) = 1.0;

%right region (x >= 0.5)
rho(x >= 0.5) = 0.125;
u(x >= 0.5) = 0.0;
p(x >= 0.5) = 0.1;

%compute initial conservative variables
E = p / (gamma - 1) + 0.5 * rho .* u.^2; %total energy
U = [rho; rho .* u; E]; %conservative variables: [density; momentum; energy]

%time-stepping loop
t = 0;
while t < t_final
    %determine time step from CFL condition
    a = sqrt(gamma * p ./ rho); %sound speed
    dt = CFL * dx / max(abs(u) + a);
    if t + dt > t_final
        dt = t_final - t; %adjust for final time step
    end

    %compute fluxes using HLLC Riemann solver
    F = zeros(size(U));
    for i = 2:Nx-1
        %left and right states
        UL = U(:, i-1);
        UR = U(:, i);

        %compute HLLC flux
        F(:, i-1) = HLLC_flux(UL, UR, gamma);
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

    %update time
    t = t + dt;
end

%plot
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

%HLLC flux function
function F = HLLC_flux(UL, UR, gamma)
    %extract left and right states
    rhoL = UL(1); uL = UL(2) / rhoL; EL = UL(3);
    rhoR = UR(1); uR = UR(2) / rhoR; ER = UR(3);

    %compute pressure and sound speed
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
