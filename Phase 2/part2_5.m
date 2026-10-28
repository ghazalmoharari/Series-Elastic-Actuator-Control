%% Solution for Phase 2: Parts 2.5 and 5
% Controller: Phase 1 PD (Kp=25, Kd=10)
clc; clear; close all;

%% --- SYSTEM PARAMETERS ---
J1 = 0.2;       % Motor Inertia
J2 = 0.5;       % Arm Inertia
k_lin = 20;     % Linear Stiffness
beta = 100;     % Non-linear Stiffness
alpha0 = 0.1;   % Operating Point
Kp = 25;        % Controller P Gain
Kd = 10;        % Controller D Gain

% Linearized Equivalent Stiffness
K_eq = k_lin + 3 * beta * alpha0^2;

%% =================================================================
%% --- PART 2.5: FREQUENCY DOMAIN ANALYSIS (Bode & Stability) ---
%% =================================================================

% 1. Define Plant Transfer Function G(s)
% G(s) = (J2*s^2 + Keq) / ( s^2 * [J1*J2*s^2 + Keq(J1+J2)] )
num_G = [J2, 0, K_eq];
den_G = [J1*J2, 0, K_eq*(J1+J2), 0, 0];
G_plant = tf(num_G, den_G);

% 2. Define Controller C(s) = Kd*s + Kp
C_ctrl = tf([Kd, Kp], [1]);

% 3. Open Loop Transfer Function L(s)
L_open = C_ctrl * G_plant;

% 4. Draw Bode Plot & Check Stability
figure('Name', 'Part 2.5: Bode Plot (Phase 1 Controller)', 'Color', 'w');
margin(L_open);
grid on;
title('Part 2.5: Bode Diagram & Stability (Kp=25, Kd=10)');

% 5. Check Closed Loop Stability
sys_cl = feedback(L_open, 1);
poles_cl = pole(sys_cl);
is_stable = all(real(poles_cl) < 0);

fprintf('--- PART 2.5 RESULTS ---\n');
if is_stable
    fprintf('System Stability: STABLE\n');
else
    fprintf('System Stability: UNSTABLE\n');
end
[Gm, Pm, Wcg, Wcp] = margin(L_open);
fprintf('Phase Margin: %.2f deg at %.2f rad/s\n', Pm, Wcp);
fprintf('Gain Margin:  %.2f dB\n\n', 20*log10(Gm));


%% =================================================================
%% --- PART 5: LOAD DISTURBANCE REJECTION ---
%% =================================================================

% 1. Simulation Settings
t_span = 0:0.001:3;
x0 = [0; 0; 0; 0]; % Start at rest

% 2. Run Non-Linear Simulation with Disturbance
[t, x] = ode45(@(t,x) nonlinear_dist_dynamics(t, x, Kp, Kd, k_lin, beta, J1, J2), t_span, x0);
theta1 = x(:,1);
theta2 = x(:,3);

% 3. Calculate Performance Metrics
% Max Deviation (from 0)
max_dev = max(abs(theta1));

% Recovery Time (Time to return to within 0.001 rad)
threshold = 0.001;
idx_settled = find(abs(theta1) > threshold, 1, 'last');
if isempty(idx_settled)
    t_rec = 0;
else
    t_rec = t(idx_settled);
end

fprintf('--- PART 5 RESULTS ---\n');
fprintf('Disturbance: Triangular (Peak 2Nm, Duration 0.5s)\n');
fprintf('Max Motor Deviation: %.4f rad\n', max_dev);
fprintf('Recovery Time:       %.4f s\n', t_rec);

% 4. Plot Results
figure('Name', 'Part 5: Disturbance Rejection', 'Color', 'w');

subplot(2,1,1);
plot(t, theta1, 'b', 'LineWidth', 1.5); hold on;
plot(t, theta2, 'r:', 'LineWidth', 1.5);
yline(threshold, 'k--', 'LineWidth', 0.5); 
yline(-threshold, 'k--', 'LineWidth', 0.5);
title('Part 5: Response to Arm Disturbance');
legend('\theta_1 (Motor)', '\theta_2 (Arm)', 'Threshold');
ylabel('Position (rad)'); grid on;

subplot(2,1,2);
% Recreate disturbance for plotting
tau_d = arrayfun(@(t) get_disturbance(t), t);
plot(t, tau_d, 'k', 'LineWidth', 1.5);
title('Applied Disturbance Torque (on Arm)');
xlabel('Time (s)'); ylabel('Torque (Nm)'); 
grid on; ylim([-0.5 2.5]);


%% =================================================================
%% --- LOCAL FUNCTIONS ---
%% =================================================================

% Non-Linear Dynamics with Disturbance
function dxdt = nonlinear_dist_dynamics(t, x, Kp, Kd, k, beta, J1, J2)
    th1 = x(1); dth1 = x(2); 
    th2 = x(3); dth2 = x(4);
    
    % Controller: Regulate to 0
    % u = -Kp*th1 - Kd*dth1
    tau_m = -Kp*th1 - Kd*dth1;
    
    % Disturbance on Arm
    tau_d = get_disturbance(t);
    
    % Spring Force
    dth = th1 - th2;
    T_spring = k*dth + beta*(dth)^3;
    
    % Equations of Motion
    % J1*ddth1 = tau_m - T_spring
    ddth1 = (tau_m - T_spring) / J1;
    
    % J2*ddth2 = T_spring + tau_d
    ddth2 = (T_spring + tau_d) / J2;
    
    dxdt = [dth1; ddth1; dth2; ddth2];
end

% Triangular Disturbance Profile
function val = get_disturbance(t)
    if t < 0
        val = 0;
    elseif t <= 0.25
        % Rise 0->2 in 0.25s (Slope 8)
        val = 8 * t;
    elseif t <= 0.5
        % Fall 2->0 in 0.25s (Slope -8)
        val = 2 - 8 * (t - 0.25);
    else
        val = 0;
    end
end