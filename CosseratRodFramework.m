function CosseratRodFramework
    %Main Function for Model Input
    
    while true
        %Display Menu
        fprintf('\nCosserat Rod Model Selection:\n');
        fprintf('1. Linear Elastic Cosserat Model with General Actuation\n');
        fprintf('2. Hyperelastic Cosserat Model with General Actuation\n');
        fprintf('3. Linear Elastic Cosserat Model with Tendon Actuation\n');
        fprintf('4. Linear Elastic Cosserat Control System with Inverse Kinematics\n');
        fprintf('5. Exit\n\n');
        
        %Model Choice
        choice = input('Select a model to run (1-5): ');
        
        %Execute Model
        switch choice
            case 1
                fprintf('\nRunning Linear Elastic Cosserat Model...\n');
                LinearElasticCosseratModel();
                
            case 2
                fprintf('\nRunning Hyperelastic Cosserat Model...\n');
                HyperelasticCosseratModel();
                
            case 3
                fprintf('\nRunning Tendon-Actuated Model...\n');
                LinearElasticTendonModel();
                
            case 4
                fprintf('\nRunning Cosserat Control System...\n');
                LinearElasticControlSystem();
                
            case 5
                fprintf('Exiting program.\n');
                return;
                
            otherwise
                fprintf('Invalid choice. Please select a number between 1 and 5.\n');
        end
    end
end

function LinearElasticCosseratModel
clc;
format shortG

%Rod/Material Properties
L = 0.4;                                    %Initial Length
E = 207e9;                                  %Young's Modulus
r = 0.0012;                                 %Radius
rho = 7850;                                 %Density
Bse = 1e-4*eye(3);                          %Damping Coefficient (Shear & Extension)
Bbt = 0.0015*eye(3);                        %Damping Coefficient (Bending and Torsion)
C = zeros(3);                               %Damping Coefficient
poisson = 0.3;                              %Poisson's Ratio

%Default Initialisation
vstar = [0; 0; 1];                          %Initial Orientation
g = 0*[-9.81; 0; 0];                        %Gravity 

%Tip Load Definition(Distal End)
F_tip = [-0.5; 0; 0];                       % Tip Force (x,y,z)
M_tip = [0; 0; 0];                          % Tip Moment (x,y,z)

%Solver Parameters
N = 200;                                    %Spatial Resolution
dt = 0.005;                                 %Time Step
Steps = 300;                                %Number of Timesteps
alpha = 0;                                  %BDF Coefficient

%Dependent Parameter Calculations
A = pi*r^2;                                 %Cross-sectional Area
J = diag([pi*r^4/4, pi*r^4/4, pi*r^4/2]);   %Inertia
G = E/(2*(1 + poisson));                    %Shear Modulus
Kse = diag([G*A, G*A, E*A]);                %Stiffness Matrix (Shear & Extension)
Kbt = diag([E*J(1,1), E*J(2,2), G*J(3,3)]); %Stiffness Matrix (Bending and Torsion)
ds = L/(N - 1);                             %Grid Spacing 

%Boundary Conditions
p0 = [0; 0; 0]; 
h0 = [1; 0; 0; 0];
q0 = [0; 0; 0]; 
w0 = [0; 0; 0];

%BDF-alpha Coefficients
c0 = (1.5 + alpha) / (dt*(1+alpha)); 
c1 = -2/dt;
c2 = (0.5 + alpha) / (dt*(1+alpha));
d1 = alpha / (1+alpha);

%Precompute Constant Terms
Kse_plus_c0_Bse_inv = inv(Kse + c0*Bse);
Kbt_plus_c0_Bbt_inv = inv(Kbt + c0*Bbt);
Kse_vstar = Kse*vstar;
rhoA = rho*A; 
rhoAg = rhoA*g; 
rhoJ = rho*J;

%Initialize Tip Displacement Storage
tip_displacement = zeros(Steps, 3); 

%Straight Configuration Initialization
y = [zeros(2, N); linspace(0, L, N); zeros(16, N)];
z = [zeros(2, N); ones(1, N); zeros(3, N)];
y_prev = y; 
z_prev = z;

%Model Solver
Dynamic_Plot;
G_guess = zeros(6, 1);                      %Shooting method initial guess


%Main Time Integration Loop
tic;
for i = 2:Steps
    %Setting History Terms
    yh = c1*y + c2*y_prev;
    zh = c1*z + c2*z_prev;
    y_prev = y;
    z_prev = z;

    %Midpoints for RK4
    yh_int = 0.5*(yh(:, 1:end-1) + yh(:, 2:end));
    zh_int = 0.5*(zh(:, 1:end-1) + zh(:, 2:end));

    %Shooting Method Call
    G_guess = fsolve(@getResidual, G_guess, optimoptions('fsolve','Display','off'));
    tip_displacement(i, :) = y(1:3, N)'; 
    Dynamic_Plot;
end
runtime = toc;

final_displacement = [y(1, :); y(2, :); y(3, :)];
disp("Final Displacements at each node:");
disp(final_displacement);
disp(['Simulation completed in ', num2str(runtime), ' seconds.'])
Tip_Displacement()

%Function Definitions
function E = getResidual(G)
    n0 = G(1:3); 
    m0 = G(4:6);
    y(:, 1) = [p0; h0; n0; m0; q0; w0];
    %Runge-Kutta Integration
    for j = 1:N-1
        yj = y(:, j); 
        yhj_int = yh_int(:, j);
        [k1, z(:, j)] = ODE(yj, yh(:, j), zh(:, j), j);
        [k2, ~] = ODE(yj + k1*ds/2, yhj_int, zh_int(:, j), j);
        [k3, ~] = ODE(yj + k2*ds/2, yhj_int, zh_int(:, j), j);
        [k4, ~] = ODE(yj + k3*ds, yh(:, j+1), zh(:, j+1), j+1);
        y(:, j + 1) = yj + ds*(k1 + 2*(k2 + k3) + k4)/6;
    end
    
    %Get Tip Boundary Conditions
    nL = y(8:10, N);
    mL = y(11:13, N);
    
    %Return Residual at Boundary
    E = [F_tip - nL; M_tip - mL];
end

function [ys, z_out] = ODE(y, yh, zh, node_idx)
    h = y(4:7); 
    n = y(8:10); 
    m = y(11:13);
    q = y(14:16); 
    w = y(17:19);
    vh = zh(1:3); 
    uh = zh(4:6);

    %Quaternion to Rotation
    h1 = h(1); 
    h2 = h(2); 
    h3 = h(3); 
    h4 = h(4);
    R = eye(3)+ 2/(h'*h)*[-h3^2 - h4^2, h2*h3 - h4*h1, h2*h4 + h3*h1;
         h2*h3 + h4*h1, -h2^2 - h4^2, h3*h4 - h2*h1;
         h2*h4 - h3*h1, h3*h4 + h2*h1, -h2^2 - h3^2];

    %Constitutive Law
    v = Kse_plus_c0_Bse_inv*(R'*n + Kse_vstar - Bse*vh);
    u = Kbt_plus_c0_Bbt_inv*(R'*m - Bbt*uh);
    z_out = [v; u];

    %Time Derivatives 
    yt = c0*y + yh;
    zt = c0*z_out + zh;
    vt = zt(1:3); 
    ut = zt(4:6);
    qt = yt(14:16); 
    wt = yt(17:19);

    %Weight and Square-Law-Drag
    f = rhoAg - R*C*q.*abs(q);
    
    %Rod State Derivatives 
    ps = R*v;
    ns = rhoA*R*(cross(w, q) + qt) - f;
    ms = R*(cross(w, rhoJ * w) + rhoJ*wt) - cross(ps, n);
    qs = vt - cross(u, q) + cross(w, v);
    ws = ut - cross(u, w);

    hs = [0, -u(1), -u(2), -u(3);
          u(1), 0, u(3), -u(2);
          u(2), -u(3), 0, u(1);
          u(3), u(2), -u(1), 0] * h / 2;

    ys = [ps; hs; ns; ms; qs; ws];
end

function Dynamic_Plot
    %Create 3D Plot
    %Reordered Coordinates (z,y,x)
    clf;
    plot3(y(3, :), y(2, :), y(1, :), 'r-', 'LineWidth', 2);  
    hold on;
    
    %Plot Tip Force if Non-Zero
    if norm(F_tip) > 0
        scale_factor = 0.2;                 
        quiver3(y(3, N), y(2, N), y(1, N), ...
               F_tip(3)*scale_factor, F_tip(2)*scale_factor, F_tip(1)*scale_factor, ...
               0, 'b', 'LineWidth', 2);
    end
    
    %Plot Tip Moment if Non-Zero
    if norm(M_tip) > 0
        scale_factor = 0.2; 
        quiver3(y(3, N), y(2, N), y(1, N), ...  
               M_tip(3)*scale_factor, M_tip(2)*scale_factor, M_tip(1)*scale_factor, ...  
               0, 'g', 'LineWidth', 2);
    end
 
    grid on;
    axis equal;
    box on;
    zlabel('x (m)');    
    ylabel('y (m)');
    xlabel('z (m)');
    zlim([-0.2*L, 0.2*L]);  
    ylim([-0.2*L, 0.2*L]);
    xlim([0, 1.1*L]);
    view(30, 30);  
    rotate3d on;
    title('Cosserat Rod Deformation');
    legend('Rod', 'Force', 'Moment');
    
    drawnow;
    pause(0.05);
end

function Tip_Displacement()
    t = (1:Steps) * dt;
    
    %X-component
    figure;
    plot(t, tip_displacement(:,1), 'r-', 'LineWidth', 1.5);
    xlabel('Time (s)');
    ylabel('x-Displacement (m)');
    title('Tip x-Displacement Over Time');
    grid on;
    
    %Y-component
    figure;
    plot(t, tip_displacement(:,2), 'g-', 'LineWidth', 1.5);
    xlabel('Time (s)');
    ylabel('y-Displacement (m)');
    title('Tip y-Displacement Over Time');
    grid on;
    
    %Z-component
    figure;
    plot(t, tip_displacement(:,3), 'b-', 'LineWidth', 1.5);
    xlabel('Time (s)');
    ylabel('z-Displacement (m)');
    title('Tip z-Displacement Over Time');
    grid on;
    
    %Magnitude
    figure;
    magnitude = sqrt(sum(tip_displacement.^2, 2));
    plot(t, magnitude, 'k-', 'LineWidth', 1.5);
    xlabel('Time (s)');
    ylabel('Displacement Magnitude (m)');
    title('Tip Displacement Magnitude');
    grid on;
end

end

function HyperelasticCosseratModel
clc;
format shortG

%Rod/Material Properties
L = 0.15;                                   %Initial Length
r = 0.003;                                  %Radius
rho = 1080;                                 %Density
Bse = 1e-5*eye(3);                          %Damping Coefficient (Shear & Extension)
Bbt = 1e-5*eye(3);                          %Damping Coefficient (Bending and Torsion)
C = zeros(3);                               %Damping Coefficient
poisson = 0.49;                             %Poisson's Ratio
mu = 0.24203*10^6;                          %Shear Modulus
K= mu/(1 - 2*poisson);                      %Bulk Modulus

%Hyperelastic Material Parameters (Neo-Hookean Model)
C10 = mu/2;                                 % Neo-Hookean parameter (μ = 2*C10)
D1 = 2/K;                                   % Compressibility parameter (K =2mu(1+v)/(3(1-2v))= 2/D1)

%Initial Modulii
E_initial = 6*C10;                          % Approximate Initial Young's Modulus (Assumption: Incompressible)
G_initial = 2*C10;                          % Initial Shear Modulus (Assumption: Incompressible)

%Default Initialisation
vstar = [0; 0; 1];                          % Initial Orientation
g = 0*[-9.81; 0; 0];                        % Gravity

%Tip Load Definition(Distal End)
F_tip = [-0.005; 0; 0];                     % Tip Force (x,y,z)
M_tip = [0; 0; 0];                          % Tip Moment (x,y,z)

%Solver Parameters
N = 200;                                    % Spatial Resolution
dt = 0.005;                                 % Time Step
Steps = 300;                                % Number of Timesteps
alpha = 0;                                  %BDF Coefficient

%Dependent Parameter Calculations
A = pi*r^2;                                 % Cross-sectional Area
J = diag([pi*r^4/4, pi*r^4/4, pi*r^4/2]);   % Inertia
                                            % Initial Stiffness (Shear & Extension)
Kse_initial = diag([G_initial*A, G_initial*A, E_initial*A]); 
                                            % Initial Stiffness (Bending and Torsion)
Kbt_initial = diag([E_initial*J(1,1), E_initial*J(2,2), G_initial*J(3,3)]);                                             
ds = L/(N - 1);                             % Grid Spacing

%Boundary Conditions 
p0 = [0; 0; 0]; 
h0 = [1; 0; 0; 0];
q0 = [0; 0; 0]; 
w0 = [0; 0; 0];

%BDF-alpha Coefficients
c0 = (1.5 + alpha) / (dt*(1+alpha)); 
c1 = -2/dt;
c2 = (0.5 + alpha) / (dt*(1+alpha));
d1 = alpha / (1+alpha);

%Precompute Constant Terms
Kse_vstar = Kse_initial*vstar;
rhoA = rho*A; 
rhoAg = rhoA*g; 
rhoJ = rho*J;

%Initialize Tip Displacement Storage
tip_displacement = zeros(Steps, 1); 

%Initialize Displacement Storage
all_displacements = zeros(N, 3, Steps);

%Initialize Strain Energy Storage
strain_energy = zeros(Steps, 1);

%Initialize Stiffness Tracking
E_values = zeros(Steps, 1);
G_values = zeros(Steps, 1);
E_values(1) = E_initial;
G_values(1) = G_initial;

%Straight Configuration Initialization
y = [zeros(2, N); linspace(0, L, N); zeros(16, N)];
z = [zeros(2, N); ones(1, N); zeros(3, N)];
y_prev = y; 
z_prev = z;

%Declaration in Global Scope for ODE Function Access
y_global = y;
z_global = z;
yh_global = zeros(size(y));
zh_global = zeros(size(z));
yh_int_global = zeros(size(yh_global, 1), size(yh_global, 2)-1);
zh_int_global = zeros(size(zh_global, 1), size(zh_global, 2)-1);

%Model Solver
Dynamic_Plot();
G_guess = zeros(6, 1);                      %Shooting Method Initial Guess

%Main Time Integration Loop
tic;
for i = 2:Steps
    % Setting History Terms
    yh_global = c1*y_global + c2*y_prev;
    zh_global = c1*z_global + c2*z_prev;
    y_prev = y_global;
    z_prev = z_global;

    % Midpoints for RK4
    yh_int_global = 0.5*(yh_global(:, 1:end-1) + yh_global(:, 2:end));
    zh_int_global = 0.5*(zh_global(:, 1:end-1) + zh_global(:, 2:end));

    % Shooting Method Call (Safety Loops)
    options = optimoptions('fsolve', 'Display', 'off', 'MaxFunctionEvaluations', 1000, ...
                          'MaxIterations', 100, 'FunctionTolerance', 1e-6);
    try
        G_guess = fsolve(@getResidual, G_guess, options);
    catch ME
        fprintf('Error at timestep %d: %s\n', i, ME.message);
        % Fall back to Previous Solution with Small Change
        G_guess = G_guess + 1e-6*randn(size(G_guess));
        % Final Try with Relaxed Options
        try
            options.FunctionTolerance = 1e-4;
            G_guess = fsolve(@getResidual, G_guess, options);
        catch
            fprintf('Still failing at timestep %d, using previous solution\n', i);
        end
    end
    
    tip_displacement(i) = y_global(1, N);
    
    % Strain Energy Calculation
    [strain_energy(i), E_avg, G_avg] = calculateStrainEnergy();
    E_values(i) = E_avg;
    G_values(i) = G_avg;
    
    % Displacement Storage
    all_displacements(:, 1, i) = y_global(1, :)'; % x-displacement
    all_displacements(:, 2, i) = y_global(2, :)'; % y-displacement
    all_displacements(:, 3, i) = y_global(3, :)'; % z-position
    Dynamic_Plot();
end
runtime = toc;

final_displacement = [y_global(1, :); y_global(2, :); y_global(3, :)]';
disp("Final Displacements at each node:");
disp(final_displacement);
disp(['Simulation completed in ', num2str(runtime), ' seconds.'])
Tip_Displacement();
Strain_Energy();
Stiffness_Variation();

%Function Definitions
function E = getResidual(G)
    n0 = G(1:3); 
    m0 = G(4:6);
    y_global(:, 1) = [p0; h0; n0; m0; q0; w0];
    
    %Runge-Kutta Integration (with Safety Loops)
    for j = 1:N-1
        yj = y_global(:, j); 
        yhj = yh_global(:, j);
        zhj = zh_global(:, j);
        yhj_int = yh_int_global(:, j);
        zhj_int = zh_int_global(:, j);     
        try
            [k1, z_global(:, j)] = ODE(yj, yhj, zhj, j);
            [k2, ~] = ODE(yj + k1*ds/2, yhj_int, zhj_int, j);
            [k3, ~] = ODE(yj + k2*ds/2, yhj_int, zhj_int, j);
            [k4, ~] = ODE(yj + k3*ds, yh_global(:, j+1), zh_global(:, j+1), j+1);
            y_global(:, j + 1) = yj + ds*(k1 + 2*(k2 + k3) + k4)/6;     
        catch ME
            %Integration Error Handling
            fprintf('Integration error at node %d: %s\n', j, ME.message);
            if j > 1
                %Forward Euler as Fallback
                try
                    [k1, z_global(:, j)] = ODE(yj, yhj, zhj, j);
                    y_global(:, j + 1) = yj + ds*k1;
                catch
                    %Extrapolate from Previous Node if Failure Persists
                    if j > 2
                        y_global(:, j + 1) = 2*y_global(:, j) - y_global(:, j-1);
                        z_global(:, j) = 2*z_global(:, j-1) - z_global(:, j-2);
                    else
                        %No Previous History, Copy Current State
                        y_global(:, j + 1) = yj;
                        z_global(:, j) = zeros(size(z_global, 1), 1);
                    end
                end
            else
                %Apply Small Change to First Node
                y_global(:, j + 1) = yj + 1e-6*randn(size(yj));
                z_global(:, j) = zeros(size(z_global, 1), 1);
            end
        end
    end
    
    %Get Tip Boundary Conditions
    nL = y_global(8:10, N);
    mL = y_global(11:13, N);
    
    %Return Residual at Boundary
    E = [F_tip - nL; M_tip - mL];
end

function [ys, z_out] = ODE(y, yh, zh, node_idx)
    h = y(4:7); 
    n = y(8:10); 
    m = y(11:13);
    q = y(14:16); 
    w = y(17:19);
    vh = zh(1:3); 
    uh = zh(4:6);

    % Check Quaternion Norm
    h_norm = norm(h);
    if h_norm < 1e-8
        h = [1; 0; 0; 0];                   % Reset to Initial Orientation if Corrupted
    elseif abs(h_norm - 1) > 1e-6
        h = h / h_norm;                     % Normalize Quaternion
    end

    %Quaternion to Rotation
    h1 = h(1); 
    h2 = h(2); 
    h3 = h(3); 
    h4 = h(4);
    R = eye(3) + 2/(h'*h)*[-h3^2 - h4^2, h2*h3 - h4*h1, h2*h4 + h3*h1;
         h2*h3 + h4*h1, -h2^2 - h4^2, h3*h4 - h2*h1;
         h2*h4 - h3*h1, h3*h4 + h2*h1, -h2^2 - h3^2];

    %Constitutive Law
    %Compute Deformation Measures
    v_current = R'*q;
    
    %Prevent Division by Zero
    if norm(vstar) < 1e-8
        J_vol = 1.0;
    else
        J_vol = max(0.5, min(2.0, norm(v_current) / norm(vstar))); % Limit extreme deformations
    end
    
    %Calculate Principal Direction Stretch
    if norm(vstar(1:2)) < 1e-8
        lambda1 = 1.0;
    else
        lambda1 = max(0.5, min(2.0, norm(v_current(1:2)) / norm(vstar(1:2)))); 
    end
    
    if abs(vstar(3)) < 1e-8
        lambda3 = 1.0;
    else
        lambda3 = max(0.5, min(2.0, v_current(3) / vstar(3)));
    end
    
    %Update Shear Modulus
    mu_current = C10 * (1 + 1/J_vol);
    
    %Update Bulk Modulus
    K_current = 2/D1 * (1 + (J_vol - 1)^2);
    
    %Compute Young's modulus from μ and K
    E_current = 9 * K_current * mu_current / (3 * K_current + mu_current);
    G_current = mu_current;
    
    %Limit Extreme Values for Stability
    E_current = max(0.5*E_initial, min(2.0*E_initial, E_current));
    G_current = max(0.5*G_initial, min(2.0*G_initial, G_current));
    
    %Update Stiffness Matrices
    Kse_current = diag([G_current*A, G_current*A, E_current*A]);
    Kbt_current = diag([E_current*J(1,1), E_current*J(2,2), G_current*J(3,3)]);
    
    % Use Blend of Initial and Current Stiffnesses for Stability
    alpha = 0.7;                            % Blend Factor (0 = use initial, 1 = use current)
    Kse_blend = alpha * Kse_current + (1-alpha) * Kse_initial;
    Kbt_blend = alpha * Kbt_current + (1-alpha) * Kbt_initial;
    Kse_plus_c0_Bse = Kse_blend + c0*Bse;
    Kbt_plus_c0_Bbt = Kbt_blend + c0*Bbt;
    
    %Safe Matrix Inversion
    try
        Kse_plus_c0_Bse_inv = inv(Kse_plus_c0_Bse);
    catch
        Kse_plus_c0_Bse_inv = inv(Kse_plus_c0_Bse + 1e-6*eye(3));
    end
    try
        Kbt_plus_c0_Bbt_inv = inv(Kbt_plus_c0_Bbt);
    catch
        Kbt_plus_c0_Bbt_inv = inv(Kbt_plus_c0_Bbt + 1e-6*eye(3));
    end
    
    v = Kse_plus_c0_Bse_inv*(R'*n + Kse_blend*vstar - Bse*vh);
    u = Kbt_plus_c0_Bbt_inv*(R'*m - Bbt*uh);
    
    %Check for NaN or Inf
    if any(isnan(v)) || any(isinf(v))
        error('NaN or Inf detected in v calculation');
    end
    if any(isnan(u)) || any(isinf(u))
        error('NaN or Inf detected in u calculation');
    end
    
    z_out = [v; u];

    %Time Derivatives 
    yt = c0*y + yh;
    zt = c0*z_out + zh;
    vt = zt(1:3); 
    ut = zt(4:6);
    qt = yt(14:16); 
    wt = yt(17:19);

    %Weight and Square-Law-Drag
    f = rhoAg - R*C*q.*abs(q);

    %Rod State Derivatives 
    ps = R*v;
    ns = rhoA*R*(cross(w, q) + qt) - f;
    ms = R*(cross(w, rhoJ * w) + rhoJ*wt) - cross(ps, n);
    qs = vt - cross(u, q) + cross(w, v);
    ws = ut - cross(u, w);

    hs = [0, -u(1), -u(2), -u(3);
          u(1), 0, u(3), -u(2);
          u(2), -u(3), 0, u(1);
          u(3), u(2), -u(1), 0] * h / 2;

    ys = [ps; hs; ns; ms; qs; ws];
    
    %Final safety check for NaN or Inf
    if any(isnan(ys)) || any(isinf(ys))
        error('NaN or Inf detected in ODE outputs');
    end
end

function [total_energy, E_avg, G_avg] = calculateStrainEnergy()
    
    total_energy = 0;
    E_sum = 0;
    G_sum = 0;
    valid_count = 0;
    
    for j = 1:N
        %Current State at Node j
        hj = y_global(4:7, j);
        qj = y_global(14:16, j);
        
        %Skip Invalid States
        if norm(hj) < 1e-8 || any(isnan(hj)) || any(isinf(hj))
            continue;
        end
        
        %Quaternion to Rotation
        h1 = hj(1); h2 = hj(2); h3 = hj(3); h4 = hj(4);
        R = eye(3) + 2/(hj'*hj)*[-h3^2 - h4^2, h2*h3 - h4*h1, h2*h4 + h3*h1;
             h2*h3 + h4*h1, -h2^2 - h4^2, h3*h4 - h2*h1;
             h2*h4 - h3*h1, h3*h4 + h2*h1, -h2^2 - h3^2];
        
        %Compute Deformation
        v_current = R'*qj;
        
        %Calculate Volumetric Jacobian J and Check Validity
        if norm(vstar) < 1e-8 || any(isnan(v_current)) || any(isinf(v_current))
            continue;
        end
        
        J_vol = norm(v_current) / norm(vstar);
        
        %Handle Extreme Values
        if J_vol < 0.1 || J_vol > 10.0
            continue;
        end
        
        % Calculate Prinicipal Direction Stretches
        lambda1 = max(0.1, min(10.0, norm(v_current(1:2)) / max(1e-8, norm(vstar(1:2)))));
        lambda3 = max(0.1, min(10.0, v_current(3) / max(1e-8, vstar(3))));
        
        %First Invariant I1 = λ1²+λ2²+λ3²
        I1 = lambda1^2 + lambda1^2 + lambda3^2; 
        
        %Neo-Hookean Strain Energy Density
        W = C10 * (I1 - 3) + 1/D1 * (J_vol - 1)^2;
        
        %Skip Node if W is Invalid
        if isnan(W) || isinf(W) || W < 0
            continue;
        end
        
        %Add to Total Energy (Strain Energy = Density*Volume)
        volume_element = A * ds;
        total_energy = total_energy + W * volume_element;
        
        %Track Material Parameters at Node
        mu_current = C10 * (1 + 1/J_vol);
        K_current = 2/D1 * (1 + (J_vol - 1)^2);
        E_current = 9 * K_current * mu_current / (3 * K_current + mu_current);
        G_current = mu_current;
        
        %Add to Averages if Reasonable
        if E_current > 0.1*E_initial && E_current < 10*E_initial && ...
           G_current > 0.1*G_initial && G_current < 10*G_initial
            E_sum = E_sum + E_current;
            G_sum = G_sum + G_current;
            valid_count = valid_count + 1;
        end
    end
    
    %Calculate Averages
    if valid_count > 0
        E_avg = E_sum / valid_count;
        G_avg = G_sum / valid_count;
    else
        E_avg = E_initial;
        G_avg = G_initial;
    end
    
    %Safety check
    if isnan(total_energy) || isinf(total_energy)
        total_energy = 0;
    end
end

function Dynamic_Plot()
    %Create 2D Plot
    figure(1);
    plot(y_global(3, :), y_global(1, :), 'r-', 'LineWidth', 2); 
    axis([0 1.1 * L -0.55 * L 0.55 * L]);
    daspect([1 1 1]);
    title('Updated Hyperelastic Cosserat Rod');
    xlabel('z (m)'); 
    ylabel('x (m)');
    grid on; 
    drawnow;
    pause(0.01);
end

function Tip_Displacement()
    figure(2);
    plot((1:Steps) * dt, tip_displacement);
    xlabel('Time (s)');
    ylabel('Tip Displacement (m)');
    title('Tip Displacement Over Time');
    grid on;
end

function Strain_Energy()
    figure(3);
    plot((1:Steps) * dt, strain_energy);
    xlabel('Time (s)');
    ylabel('Strain Energy (J)');
    title('Total Strain Energy Over Time');
    grid on;
end

function Stiffness_Variation()
    figure(4);
    plot((1:Steps) * dt, E_values/E_initial, 'b-', 'LineWidth', 1.5);
    hold on;
    plot((1:Steps) * dt, G_values/G_initial, 'r--', 'LineWidth', 1.5);
    hold off;
    xlabel('Time (s)');
    ylabel('Normalized Stiffness');
    title('Stiffness Variation Over Time');
    legend('E/E_0', 'G/G_0');
    grid on;
end

end

function LinearElasticTendonModel
clc;
format shortG

%Rod/Material Properties
L = 0.4;                                    %Initial Length
E = 207e9;                                  %Young's Modulus
r = 0.001;                                  %Radius
rho = 8000;                                 %Density 
Bse = zeros(3);                             %Damping Coefficient (Shear & Extension)
Bbt = diag([1e-6, 1e-6, 1e-6]);             %Damping Coefficient (Bending and Torsion)
C = diag([0.03, 0.03, 0.03]);               %Damping Coefficient
poisson = 0.3;                              %Poisson's Ratio

%Default Initialisation
g = [0;0;-9.81];                            %Gravity 

%Tendon Actuation Properties
Num_Tendons = 4;                            % Number of Tendons
tau = [5; 2; 0; 0];                         % Tendon Tensions
compliance = [1e-4; 1e-4; 1e-4; 1e-4];      % Tendon Compliance
Num_Disks = 1;
Tendon_Offset = 0.01;                       % Tendon Offset

%Four Tendons at 90 degree Intervals
rot1 = @(theta) Tendon_Offset*[cos(theta); sin(theta); 0];
rpos = {rot1(0), rot1(pi/2), rot1(pi), rot1(3*pi/2)}; 

%Solver Parameters
N = 50;                                     %Spatial Resolution
dt = 0.05;                                  %Time Step
Steps = 100;                                %Numver of Timesteps
alpha = 0;                                  %BDF Coefficient
global t;                                   %Global Time Parameter
t = 0;

%Dependent Parameter Calculations
Area = pi*r^2;                              %Cross-sectional Area
J = diag([pi*r^4/4  pi*r^4/4  pi*r^4/2]);   %Inertia
G_Modulus = E/( 2*(1+poisson) );            %Shear Modulus
                                            %Stiffness Matrix (Shear & Extension)
Kse = diag([G_Modulus*Area, G_Modulus*Area, E*Area]); 
                                            %Stiffness Matrix (Bending and Torsion)
Kbt = diag([E*J(1,1), E*J(2,2), G_Modulus*J(3,3)]); 
ds = L/(N-1);                               %Grid Spacing

%Boundary Conditions 
p0 = [0;0;0];
R0 = eye(3);
q0 = [0;0;0];
w0 = [0;0;0];

%Step Input
    function displace = Z_t(t)
        displace = 0; 
    end

%BDF-alpha Coefficients
c0 = (1.5 + alpha) / (dt*(1+alpha)); 
c1 = -2/dt;
c2 = (0.5 + alpha) / (dt*(1+alpha));
d1 = alpha / (1+alpha);

%Precompute Constant Terms
rhoA = rho*Area;       
rhoAg = rho*Area*g;

%Global Parameter Declaration
global Y Z Y_prev Z_prev;

%Initialize Displacement Storage
tip_displacement = zeros(Steps, 3); 
centerline_history = zeros(3, N, Steps);  
tendon_positions_history = zeros(3*Num_Tendons, N, Steps); 

%Configuration Initialization
Y = zeros(24+Num_Tendons, N);
Y(3, :) = linspace(0, L, N); 
for i=1:N
    Y(4:12, i) = reshape(R0, 9, 1); 
end
Z = zeros(18, N);
Z(3, :) = ones(1, N);
Y_prev = Y;
Z_prev = Z;
global i_Tstep;
global Z_h;

%Model Solver
Dynamic_Plot();
G_guess = fsolve(@staticBVP, zeros(6+Num_Tendons, 1)); 

%Main Time Integration Loop
tic;
for tStep = 1 : Steps
    %Setting History Terms
    Z_h = c1*Z+c2*Z_prev;
    Y_prev = Y;
    Z_prev = Z;
    
    %Midpoints for RK4
    Zh_int = 0.5*(Z_h(:,1:end-1) + Z_h(:,2:end));

    %Shooting Method Call
    G_guess = fsolve(@dynamicBVP, G_guess);     
    Dynamic_Plot(); 
    tip_displacement(tStep,:) = Y(1:3,end)';
    centerline_history(:,:,tStep) = Y(1:3,:);
    for i = 1:N
        p_current = Y(1:3,i);
        R_current = reshape(Y(4:12,i),3,3);
        for j = 1:Num_Tendons
            tendon_positions_history(3*j-2:3*j, i, tStep) = p_current + R_current*rpos{j};
        end
    end
    
    disp(['Time: ', num2str(t), ' Tip position: ', num2str(Y(1:3,end)')]);
    t = t+dt;                               % Update Time 
end
runtime = toc;


disp(['Simulation completed in ', num2str(runtime), ' seconds.'])

%Function Definitions
function [ys] = staticODE(s,y)
    R = reshape(y(4:12),3,3);
    v = y(13:15);
    u = y(16:18);

    %Setup Linear Tendon System
    a = zeros(3,1);
    b = zeros(3,1);
    A = zeros(3,3);
    G = zeros(3,3);
    H = zeros(3,3);
    pib_s_norm = y(24+1:24+Num_Tendons);

    for i = 1 : Num_Tendons
        pb_si = cross(u,rpos{i}) + v;
        pib_s_norm(i) = norm(pb_si);
        A_i = -hat(pb_si)^2*(tau(i)/pib_s_norm(i)^3);
        G_i = -A_i*hat(rpos{i});
        a_i = A_i*cross(u,pb_si);

        a = a + a_i;
        b = b + cross(rpos{i}, a_i);
        A = A + A_i;
        G = G + G_i;
        H = H + hat(rpos{i})*G_i;
    end

    K = [A + Kse, G ; G.', H + Kbt];
    nb = Kse*(v - [0;0;1]);
    mb = Kbt*u;

    rhs = [-cross(u,nb) - R'*rho*Area*g - a; -cross(u,mb) - cross(v,nb) - b];

    %Rod State Derivatives 
    ps = R*v;
    Rs = R*hat(u);
    vs_and_us = K\rhs;
    qs = zeros(3,1);
    ws = zeros(3,1);
    ys = [ps; reshape(Rs,9,1); vs_and_us; qs; ws; pib_s_norm];
end

function [ys,z] = dynamicODE(y,zh)
    R = reshape(y(4:12),3,3);
    v = y(13:15);
    u = y(16:18);
    q = y(19:21);
    w = y(22:24);

    vh = zh(1:3);
    uh = zh(4:6);
    qh = zh(7:9);
    wh = zh(10:12);
    vsh = zh(13:15);
    ush = zh(16:18);
    
    %Setup Linear Tendon System
    a = zeros(3,1);
    b = zeros(3,1);
    A = zeros(3,3);
    G = zeros(3,3);
    H = zeros(3,3);
    pib_s_norm = zeros(Num_Tendons,1);

    for num = 1 : Num_Tendons
        pb_si = cross(u,rpos{num}) + v;
        pib_s_norm(num) = norm(pb_si);
        A_i = -hat(pb_si)^2 * (tau(num)/pib_s_norm(num)^3);
        G_i = -A_i * hat(rpos{num});
        a_i = A_i * cross(u,pb_si);

        a = a + a_i;
        b = b + cross(rpos{num}, a_i);
        A = A + A_i;
        G = G + G_i;
        H = H + hat(rpos{num})*G_i;
    end

    K = [A + Kse + c0*Bse,G ; G.',  H + Kbt + c0*Bbt];

    v_t = c0*v + vh;
    u_t = c0*u + uh;
    q_t = c0*q + qh;
    w_t = c0*w + wh;

    nb = Kse*(v - [0;0;1])+Bse*v_t;
    mb = Kbt*u+Bbt*u_t;

    rhs = [-a + rhoA*(cross(w,q)+q_t) + C*q.*norm(q) - R'*rhoAg - cross(u,nb)-Bse*vsh;
          -b + rho*cross(w,J*w)+rho*J*w_t - cross(v,nb)-cross(u,mb)-Bbt*ush];

    % ODEs
    ps = R*v;
    Rs = R*hat(u);
    vs_and_us = K\rhs;
    qs = v_t - hat(u)*q + hat(w)*v;
    ws = u_t - hat(u)*w;

    % Pack state vector derivative
    z = [v;u;q;w;vs_and_us];
    ys = [ps; reshape(Rs,9,1); vs_and_us; qs; ws; pib_s_norm];
end

function distal_error = staticBVP(guess)
    n0 = guess(1:3);
    v0 = Kse\n0 + [0;0;1];
    u0 = guess(4:6);
    tau_guess = max(guess(7:7+Num_Tendons-1),0); 
    slack = -min(guess(7:7+Num_Tendons-1),0);

    %Step Input
    pi_b_norm = -Z_t(t)*ones(Num_Tendons,1);

    s = linspace(0,L,N);

    %Euler's Method
    Y(:,1) = [p0; reshape(R0,9,1); v0; u0; q0; w0; pi_b_norm]; 
    for j = 1 : N-1
        [ys] = staticODE(s(j),Y(:,j));
        Y(:,j+1) = Y(:,j)+ds*ys;
        Z(1:6,j) = Y(13:18,j);
        Z(13:18,j) = ys(13:18);
    end

    %Internal Forces in Backbone
    vL = Y(13:15,end);
    uL = Y(16:18,end);
    nb_L = Kse*(vL - [0;0;1]);
    mb_L = Kbt*uL;

    %Equilibrium Error at Tip
    force_error = -nb_L;
    moment_error = -mb_L;

    for index = 1 : Num_Tendons
        pb_si = cross(uL,rpos{index}) + vL;
        Fb_i = -tau(index)*pb_si/norm(pb_si);
        force_error = force_error + Fb_i;
        moment_error = moment_error + cross(rpos{index}, Fb_i);
    end

    %Length Violation Error
    integrated_lengths = Y(24+1:24+Num_Tendons,end);
    l_star = L;
    stretch = l_star*(compliance.*tau_guess); % calculate first for fast
    length_error = (integrated_lengths + slack)-(l_star + stretch);
    distal_error = [force_error; moment_error; length_error];
end

function distal_error = dynamicBVP(guess)
    n0 = guess(1:3);
    v0 = Kse\n0 + [0;0;1];
    u0 = guess(4:6);
    tau_guess = max(guess(7:7+Num_Tendons-1),0); 
    slack = -min(guess(7:7+Num_Tendons-1),0);

    %Step Input
    pi_b_norm = -Z_t(t)*ones(Num_Tendons,1);

    %Euler's method
    Y(:,1) = [p0; reshape(R0,9,1); v0; u0; q0; w0; pi_b_norm]; 
    for j = 1 : N-1
        [ys,Z(:,j)] = dynamicODE(Y(:,j),Z_h(:,j));
        Y(:,j+1) = Y(:,j)+ds*ys;
    end

    
    vL = Y(13:15,end);
    uL = Y(16:18,end);
    vL_t = c0*vL + Z_h(1:3,end);
    uL_t = c0*uL + Z_h(4:6,end);

    nb_L = Kse*(vL - [0;0;1])+Bse*vL_t;
    mb_L = Kbt*uL + Bbt*uL_t;

    %Internal Forces in Backbone
    force_error = -nb_L;
    moment_error = -mb_L;

    %Equilibrium Error 
    for index = 1 : Num_Tendons
        pb_si = cross(uL,rpos{index}) + vL;
        Fb_i = -tau(index)*pb_si/norm(pb_si);
        force_error = force_error + Fb_i;
        moment_error = moment_error + cross(rpos{index}, Fb_i);
    end

    %Length Violation Error
    integrated_lengths = Y(24+1:24+Num_Tendons,end);
    l_star = L;
    stretch = l_star*(compliance.*tau_guess); % calculate first for fast
    length_error = integrated_lengths + slack-(l_star + stretch);
    distal_error = [force_error; moment_error; length_error];
end

function skew_symmetric_matrix = hat(y)
    skew_symmetric_matrix = [  0   -y(3)  y(2) ;
                             y(3)   0   -y(1) ;
                            -y(2)  y(1)   0  ];
end

function Dynamic_Plot()
    centerline = Y(1:3,:);
    tendonlines = zeros(3*Num_Tendons, size(centerline,2));
    for i = 1 : size(centerline,2)
        p_show = Y(1:3,i);
        R_show = reshape(Y(4:12,i),3,3);
        for j = 1 : Num_Tendons
            tendonlines(3*j-2 : 3*j, i) = p_show + R_show*rpos{j};
        end
    end
    disks = zeros(3,4*Num_Disks);
    for i = 1 : Num_Disks
        j = round(size(centerline,2) * i / Num_Disks);
        p_show = Y(1:3,j);
        R_show = reshape(Y(4:12,j),3,3);
        disks(1:3, 4*i-3:4*i-1) = R_show;
        disks(1:3, 4*i) = p_show;
    end
    Visualize(centerline,tendonlines,Num_Tendons, disks,Num_Disks);
end

function Visualize(centerline, tendonlines, Num_Tendons, disks, Num_Disks)
    figure(1); 
    clf;
    hold on;

    %Plot Centerline
    plot3(centerline(1,:), centerline(2,:), centerline(3,:), 'k-', 'LineWidth', 2);

    %Plot Tendons
    tendon_colors = {'r-', 'g-', 'b-', 'c-'};
    for i = 1:Num_Tendons
        tendon_i = tendonlines(3*i-2:3*i, :);
        plot3(tendon_i(1,:), tendon_i(2,:), tendon_i(3,:), tendon_colors{mod(i-1,4)+1}, 'LineWidth', 1);
    end

    %Plot Disks
    for i = 1:Num_Disks
        disk_i = disks(:, 4*i-3:4*i);
        R = disk_i(:, 1:3);
        p = disk_i(:, 4);

        %Disc Visualization
        theta = linspace(0, 2*pi, 30);
        r_disk = 0.004; 
        disk_x = r_disk * cos(theta);
        disk_y = r_disk * sin(theta);
        disk_z = zeros(size(theta));
        disk_points = zeros(3, length(theta));
        for j = 1:length(theta)
            point = [disk_x(j); disk_y(j); disk_z(j)];
            disk_points(:, j) = p + R * point;
        end

        %Plot disk
        fill3(disk_points(1,:), disk_points(2,:), disk_points(3,:), 'y', 'FaceAlpha', 0.5);
    end

    grid on;
    axis equal;
    xlabel('x (m)');
    ylabel('y (m)');
    zlabel('z (m)');
    title('Tendon-Driven Cosserat Rod Deformation');
    view(3);
    drawnow;
end

end

function LinearElasticControlSystem
clc;
format shortG

%Target Displacement(x,y,z)
target_displacement = [-0.0063; 0.004; 0.39994]; 

%Controller Parameters
max_iterations = 10;                        %Maximum Iterations for Force Adjustment
convergence_tolerance = 1e-3;               %Convergence Criterion for Displacement Error
steady_state_tolerance = 1e-5;              %Tolerance for Steady-state Detection

%PID Controller Gains
K_p = 0.7;                                  %Proportional Gain
K_i = 0.5;                                  %Integral Gain
K_d = 0;                                    %Derivative Gain

%Initial Guess(x,y,z)
F_tip = [-0.05; -0.05; -0.001];          

%History Storage
[final_force, final_displacement, error_history, force_history, pid_history] = runPIDController(target_displacement, F_tip, max_iterations, convergence_tolerance, steady_state_tolerance, K_p, K_i, K_d);

disp('Inverse Kinematics Results:');
disp('-------------------------');
disp(['Target displacement: [', num2str(target_displacement(1)), ', ', num2str(target_displacement(2)), ', ', num2str(target_displacement(3)), ']']);
disp(['Final displacement:  [', num2str(final_displacement(1)), ', ', num2str(final_displacement(2)), ', ', num2str(final_displacement(3)), ']']);
disp(['Final force:         [', num2str(final_force(1)), ', ', num2str(final_force(2)), ', ', num2str(final_force(3)), ']']);
disp(['Final error:         ', num2str(norm(final_displacement - target_displacement))]);
disp(['Required iterations: ', num2str(length(error_history))]);
plotEnhancedResults(error_history, force_history, target_displacement, pid_history);

function [final_force, final_displacement, error_history, force_history, pid_history] = runPIDController(target_displacement, initial_force, max_iterations, convergence_tolerance, steady_state_tolerance, K_p, K_i, K_d)
    %Initialize Storage
    error_history = zeros(max_iterations, 1);
    force_history = zeros(max_iterations, 3);
    pid_history = struct();
    pid_history.P_terms = zeros(max_iterations, 3);
    pid_history.I_terms = zeros(max_iterations, 3);
    pid_history.D_terms = zeros(max_iterations, 3);
    pid_history.total_terms = zeros(max_iterations, 3);
    pid_history.displacements = zeros(max_iterations, 3);
    
    %Initial Force
    F_tip = initial_force;
    
    %Initialize PID Controller Variables
    error_integral = zeros(3, 1);           %Integral Term
    previous_error = zeros(3, 1);           %For Derivative Term
    previous_displacement = [];             %Store Previous Displacement for Derivative Term
    
    %Main Loop
    for iter = 1:max_iterations
        disp(['Iteration ', num2str(iter), '/', num2str(max_iterations)]);
        
        % Run Forward Model with Current Forces
        [current_displacement, convergence_iterations] = runForwardModel(F_tip, steady_state_tolerance);
        
        %Store Displacement History
        pid_history.displacements(iter, :) = current_displacement';
        
        %Error Calculation
        displacement_error = target_displacement - current_displacement;
        error_norm = norm(displacement_error);
        error_history(iter) = error_norm;
        force_history(iter, :) = F_tip';
        
        % Display current status
        disp(['  Force: [', num2str(F_tip(1)), ', ', num2str(F_tip(2)), ', ', num2str(F_tip(3)), ']']);
        disp(['  Displacement: [', num2str(current_displacement(1)), ', ', num2str(current_displacement(2)), ', ', num2str(current_displacement(3)), ']']);
        disp(['  Error: ', num2str(error_norm)]);
        disp(['  Forward model converged in ', num2str(convergence_iterations), ' steps']);
        
        %Tolerance Check
        if error_norm < convergence_tolerance
            disp('Target displacement achieved within tolerance!');
            break;
        end
        
        %PID Controller Implementation
        
        %Proportional Term
        P_term = K_p * displacement_error;
        pid_history.P_terms(iter, :) = P_term';
        
        %Integral Term (with Anti-Windup)
        %Only Integrate Controller is not Saturated
        if iter > 1 && error_history(iter) < error_history(iter-1)
            error_integral = error_integral + displacement_error;
        else
            %Reset Integral Term if Error is Increasing
            error_integral = 0.5 * error_integral;
        end
        I_term = K_i * error_integral;
        pid_history.I_terms(iter, :) = I_term';
        
        %Derivative term 
        if ~isempty(previous_displacement)
            D_term = -K_d * (current_displacement - previous_displacement);
        else
            D_term = zeros(3, 1);
        end
        previous_displacement = current_displacement;
        previous_error = displacement_error;  
        pid_history.D_terms(iter, :) = D_term';
        
        %Compute Total Control Action
        force_update = P_term + I_term + D_term;
        pid_history.total_terms(iter, :) = force_update';
        
        %Limit Maximum Force Change per Iteration
        max_force_change = 0.2;  
        for i = 1:3
            if abs(force_update(i)) > max_force_change
                force_update(i) = sign(force_update(i)) * max_force_change;
            end
        end
        
        %Update Force
        F_tip = F_tip + force_update;
        
        %Damping Factor to Prevent Oscillations
        if iter > 2 && error_history(iter) > error_history(iter-1) && error_history(iter-1) > error_history(iter-2)
            F_tip = 0.7 * F_tip + 0.3 * force_history(iter-1, :)';
            disp('  Applied damping to force update due to increasing error');
        end
        
        disp(['  P-term: [', num2str(P_term(1)), ', ', num2str(P_term(2)), ', ', num2str(P_term(3)), ']']);
        disp(['  I-term: [', num2str(I_term(1)), ', ', num2str(I_term(2)), ', ', num2str(I_term(3)), ']']);
        disp(['  D-term: [', num2str(D_term(1)), ', ', num2str(D_term(2)), ', ', num2str(D_term(3)), ']']);
    end
    
    %Return Final Values
    final_force = F_tip;
    final_displacement = current_displacement;
    
    %Trim Unused Entries in Arrays
    error_history = error_history(1:iter);
    force_history = force_history(1:iter, :);
    
    %Trim PID History
    pid_history.P_terms = pid_history.P_terms(1:iter, :);
    pid_history.I_terms = pid_history.I_terms(1:iter, :);
    pid_history.D_terms = pid_history.D_terms(1:iter, :);
    pid_history.total_terms = pid_history.total_terms(1:iter, :);
    pid_history.displacements = pid_history.displacements(1:iter, :);
end

function [final_displacement, iterations] = runForwardModel(F_tip, steady_state_tolerance)
    %Rod/Material Properties
    L = 0.4;                                %Initial Length
    E = 207e9;                              %Young's Modulus
    r = 0.0012;                             %Radius
    rho = 7850;                             %Density
    Bse = 1e-4*eye(3);                      %Damping Coefficient (Shear & Extension)
    Bbt = 0.0015*eye(3);                    %Damping Coefficient (Bending and Torsion)
    C = zeros(3);                           %Damping Coefficient
    poisson = 0.3;                          %Poisson's Ratio

    %Default Initialisation
    vstar = [0; 0; 1];                      %Initial Orientation
    g = 0*[-9.81; 0; 0];                    %Gravity 
    
    %Tip Load Definition(Distal End)
    M_tip = [0; 0; 0];                      %Moment applied at tip (x,y,z)

    %Solver Parameters
    N = 200;                                %Spatial Resolution
    dt = 0.005;                             %Time Step
    Steps = 600;                            %Number of Timesteps
    alpha = 0;                              %BDF Coefficient
    %Dependent Parameter Calculations
    A = pi*r^2;                             %Cross-sectional Area
                                            %Inertia
    J = diag([pi*r^4/4, pi*r^4/4, pi*r^4/2]); 
    G = E/(2*(1 + poisson));                %Shear Modulus
    Kse = diag([G*A, G*A, E*A]);            %Stiffness Matrix (Shear & Extension)
                                            %Stiffness Matrix (Bending and Torsion)
    Kbt = diag([E*J(1,1), E*J(2,2), G*J(3,3)]); 
    ds = L/(N - 1);                         %Grid Spacing 

    %Boundary Conditions
    p0 = [0; 0; 0]; 
    h0 = [1; 0; 0; 0];
    q0 = [0; 0; 0]; 
    w0 = [0; 0; 0];
 
    %BDF-alpha Coefficients
    c0 = (1.5 + alpha) / (dt*(1+alpha)); 
    c1 = -2/dt;
    c2 = (0.5 + alpha) / (dt*(1+alpha));
    d1 = alpha / (1+alpha);

    %Precompute Constant Terms
    Kse_plus_c0_Bse_inv = inv(Kse + c0*Bse);
    Kbt_plus_c0_Bbt_inv = inv(Kbt + c0*Bbt);
    Kse_vstar = Kse*vstar;
    rhoA = rho*A; 
    rhoAg = rhoA*g; 
    rhoJ = rho*J;

    %Initialize Tip Displacement Storage
    tip_displacement_history = zeros(Steps, 3);
    
    %Straight Configuration Initialization
    y = [zeros(2, N); linspace(0, L, N); zeros(16, N)];
    z = [zeros(2, N); ones(1, N); zeros(3, N)];
    y_prev = y; 
    z_prev = z;

    %Model Solver
    G_guess = zeros(6, 1);                  %Shooting Method Initial Guess

    %Main Time Integration Loop
    for i = 2:Steps
        %Setting History Terms
        yh = c1*y + c2*y_prev;
        zh = c1*z + c2*z_prev;
        y_prev = y;
        z_prev = z;

        %Midpoints for RK4
        yh_int = 0.5*(yh(:, 1:end-1) + yh(:, 2:end));
        zh_int = 0.5*(zh(:, 1:end-1) + zh(:, 2:end));

        %Shooting Method Call
        G_guess = fsolve(@getResidual, G_guess, optimoptions('fsolve','Display','off'));
        tip_displacement_history(i, :) = y(1:3, N)'; 

        % Check for steady state (convergence)
        if i > 10
            recent_displacements = tip_displacement_history(i-10:i, :);
            displacement_variation = max(max(abs(diff(recent_displacements))));
            if displacement_variation < steady_state_tolerance
                iterations = i;
                break;
            end
        end
    end
    if exist('iterations', 'var') == 0
        iterations = Steps;
    end
    
    
    final_displacement = y(1:3, N);
    
    %Nested Function Definiton
    function E = getResidual(G)
        n0 = G(1:3); 
        m0 = G(4:6);
        y(:, 1) = [p0; h0; n0; m0; q0; w0];
        %Runge-Kutta Integration
        for j = 1:N-1
            yj = y(:, j); 
            yhj_int = yh_int(:, j);
            [k1, z(:, j)] = ODE(yj, yh(:, j), zh(:, j), j);
            [k2, ~] = ODE(yj + k1*ds/2, yhj_int, zh_int(:, j), j);
            [k3, ~] = ODE(yj + k2*ds/2, yhj_int, zh_int(:, j), j);
            [k4, ~] = ODE(yj + k3*ds, yh(:, j+1), zh(:, j+1), j+1);
            y(:, j + 1) = yj + ds*(k1 + 2*(k2 + k3) + k4)/6;
        end
        
        %Get Tip Boundary Conditions
        nL = y(8:10, N);
        mL = y(11:13, N);
        
        %Return Residual at Boundary
        E = [F_tip - nL; M_tip - mL];
    end

    function [ys, z_out] = ODE(y, yh, zh, node_idx)
        h = y(4:7); 
        n = y(8:10); 
        m = y(11:13);
        q = y(14:16); 
        w = y(17:19);
        vh = zh(1:3); 
        uh = zh(4:6);

        %Quaternion to Rotation
        h1 = h(1); 
        h2 = h(2); 
        h3 = h(3); 
        h4 = h(4);
        R = eye(3)+ 2/(h'*h)*[-h3^2 - h4^2, h2*h3 - h4*h1, h2*h4 + h3*h1;
             h2*h3 + h4*h1, -h2^2 - h4^2, h3*h4 - h2*h1;
             h2*h4 - h3*h1, h3*h4 + h2*h1, -h2^2 - h3^2];

        %Constitutive Law
        v = Kse_plus_c0_Bse_inv*(R'*n + Kse_vstar - Bse*vh);
        u = Kbt_plus_c0_Bbt_inv*(R'*m - Bbt*uh);
        z_out = [v; u];

        %Time Derivatives 
        yt = c0*y + yh;
        zt = c0*z_out + zh;
        vt = zt(1:3); 
        ut = zt(4:6);
        qt = yt(14:16); 
        wt = yt(17:19);

        %Weight and Square-Law-Drag
        f = rhoAg - R*C*q.*abs(q);
        
        %Rod State Derivatives 
        ps = R*v;
        ns = rhoA*R*(cross(w, q) + qt) - f;
        ms = R*(cross(w, rhoJ * w) + rhoJ*wt) - cross(ps, n);
        qs = vt - cross(u, q) + cross(w, v);
        ws = ut - cross(u, w);

        hs = [0, -u(1), -u(2), -u(3);
              u(1), 0, u(3), -u(2);
              u(2), -u(3), 0, u(1);
              u(3), u(2), -u(1), 0] * h / 2;

        ys = [ps; hs; ns; ms; qs; ws];
    end
end

function plotEnhancedResults(error_history, force_history, target_displacement, pid_history)
    iteration_indices = 1:length(error_history);
    
    %Error Convergence Plot
    figure('Name', 'Error Convergence', 'Position', [100, 100, 800, 600]);
    subplot(2,1,1);
    semilogy(iteration_indices, error_history, 'b-o', 'LineWidth', 1.5);
    grid on;
    xlabel('Iteration');
    ylabel('Displacement Error (L2 norm)');
    title('Convergence of Tip Displacement Error (Log Scale)');
    
    subplot(2,1,2);
    plot(iteration_indices, error_history, 'b-o', 'LineWidth', 1.5);
    grid on;
    xlabel('Iteration');
    ylabel('Displacement Error (L2 norm)');
    title('Convergence of Tip Displacement Error (Linear Scale)');
    
    %Force Evolution Plot
    figure('Name', 'Force Evolution', 'Position', [150, 150, 800, 600]);
    plot(iteration_indices, force_history(:,1), 'r-o', 'LineWidth', 1.5);
    hold on;
    plot(iteration_indices, force_history(:,2), 'g-o', 'LineWidth', 1.5);
    plot(iteration_indices, force_history(:,3), 'b-o', 'LineWidth', 1.5);
    grid on;
    xlabel('Iteration');
    ylabel('Force Component (N)');
    title('Evolution of Force Components');
    legend('F_x', 'F_y', 'F_z');
    
    %PID Component Plot
    figure('Name', 'PID Components', 'Position', [200, 200, 1000, 800]);
    
    %X direction
    subplot(3,1,1);
    plot(iteration_indices, pid_history.P_terms(:,1), 'r-', 'LineWidth', 1.5);
    hold on;
    plot(iteration_indices, pid_history.I_terms(:,1), 'g-', 'LineWidth', 1.5);
    plot(iteration_indices, pid_history.D_terms(:,1), 'b-', 'LineWidth', 1.5);
    plot(iteration_indices, pid_history.total_terms(:,1), 'k--', 'LineWidth', 2);
    grid on;
    title('PID Components - x Direction');
    xlabel('Iteration');
    ylabel('Control Action');
    legend('P Term', 'I Term', 'D Term', 'Total');
    
    %Y direction
    subplot(3,1,2);
    plot(iteration_indices, pid_history.P_terms(:,2), 'r-', 'LineWidth', 1.5);
    hold on;
    plot(iteration_indices, pid_history.I_terms(:,2), 'g-', 'LineWidth', 1.5);
    plot(iteration_indices, pid_history.D_terms(:,2), 'b-', 'LineWidth', 1.5);
    plot(iteration_indices, pid_history.total_terms(:,2), 'k--', 'LineWidth', 2);
    grid on;
    title('PID Components - y Direction');
    xlabel('Iteration');
    ylabel('Control Action');
    legend('P Term', 'I Term', 'D Term', 'Total');
    
    %Z direction
    subplot(3,1,3);
    plot(iteration_indices, pid_history.P_terms(:,3), 'r-', 'LineWidth', 1.5);
    hold on;
    plot(iteration_indices, pid_history.I_terms(:,3), 'g-', 'LineWidth', 1.5);
    plot(iteration_indices, pid_history.D_terms(:,3), 'b-', 'LineWidth', 1.5);
    plot(iteration_indices, pid_history.total_terms(:,3), 'k--', 'LineWidth', 2);
    grid on;
    title('PID Components - z Direction');
    xlabel('Iteration');
    ylabel('Control Action');
    legend('P Term', 'I Term', 'D Term', 'Total');
    
    %Displacement Trajectory
    figure('Name', 'Displacement Trajectory', 'Position', [250, 250, 900, 700]);
    
    %X displacement
    subplot(3,1,1);
    plot(iteration_indices, pid_history.displacements(:,1), 'b-o', 'LineWidth', 1.5);
    hold on;
    plot([1, length(error_history)], [target_displacement(1), target_displacement(1)], 'r--', 'LineWidth', 1.5);
    grid on;
    xlabel('Iteration');
    ylabel('x Displacement');
    title('x Displacement vs. Target');
    legend('Actual', 'Target');
    
    %Y displacement
    subplot(3,1,2);
    plot(iteration_indices, pid_history.displacements(:,2), 'b-o', 'LineWidth', 1.5);
    hold on;
    plot([1, length(error_history)], [target_displacement(2), target_displacement(2)], 'r--', 'LineWidth', 1.5);
    grid on;
    xlabel('Iteration');
    ylabel('y Displacement');
    title('y Displacement vs. Target');
    legend('Actual', 'Target');
    
    %Z displacement
    subplot(3,1,3);
    plot(iteration_indices, pid_history.displacements(:,3), 'b-o', 'LineWidth', 1.5);
    hold on;
    plot([1, length(error_history)], [target_displacement(3), target_displacement(3)], 'r--', 'LineWidth', 1.5);
    grid on;
    xlabel('Iteration');
    ylabel('z Displacement');
    title('z Displacement vs. Target');
    legend('Actual', 'Target');
end

end
