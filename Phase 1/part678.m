%% SEA Project: Parts 6, 7, and 8
% Controller: PD (Kp=25, Kd=10)
clc; clear; close all;

%% --- SYSTEM PARAMETERS ---
J1 = 0.2;       % Motor Inertia
J2 = 0.5;       % Arm Inertia
k_lin = 20;     % Linear Stiffness
beta = 100;     % Non-linear Stiffness
alpha0 = 0.1;   % Operating Point
Kp = 25;        % Proportional Gain
Kd = 10;        % Derivative Gain

% Linearized Equivalent Stiffness
K_eq = k_lin + 3 * beta * alpha0^2;

%% --- PART 6: LINEAR MODEL STEP RESPONSE ---
% Transfer Function G(s)
num = [J2, 0, K_eq];
den = [J1*J2, 0, K_eq*(J1+J2), 0, 0];
G_plant = tf(num, den);

% Closed Loop System
C_ctrl = pid(Kp, 0, Kd);
sys_cl = feedback(C_ctrl * G_plant, 1);

% Metrics
info_lin = stepinfo(sys_cl);

% FIGURE 1: PART 6
figure('Name', 'Part 6: Linear Model Step Response', 'Color', 'w');
step(sys_cl);
grid on;
title(sprintf('Part 6: Linear Step Response (Kp=%.0f, Kd=%.0f)', Kp, Kd));
ylabel('Position (rad)');

% Output to Command Window
fprintf('--- PART 6: LINEAR DESIGN ---\n');
fprintf('Settling Time: %.4f s\n', info_lin.SettlingTime);
fprintf('Overshoot:     %.2f %%\n\n', info_lin.Overshoot);

%% --- PART 7: LINEAR VS NON-LINEAR COMPARISON ---
t_span = 0:0.01:6;
x0 = [0; 0; 0; 0]; % Rest conditions

% 1. Linear Simulation (using TF step command data)
[y_lin, t_lin] = step(sys_cl, t_span);

% 2. Non-Linear Simulation (ODE45)
[t_nl, x_nl] = ode45(@(t,x) nonlinear_dynamics(t, x, Kp, Kd, k_lin, beta, J1, J2), t_span, x0);
theta1_nl = x_nl(:,1);

% Calculate Non-Linear Metrics
final_val = 1.0;
tol = 0.02 * final_val;
is_settled = abs(theta1_nl - final_val) < tol;
idx_unsettled = find(~is_settled, 1, 'last');
if isempty(idx_unsettled), Ts_nl = 0; else, Ts_nl = t_nl(idx_unsettled); end
OS_nl = (max(theta1_nl) - final_val)/final_val * 100;

% FIGURE 2: PART 7
figure('Name', 'Part 7: Linear vs Non-Linear', 'Color', 'w');
plot(t_lin, y_lin, 'b--', 'LineWidth', 1.5); hold on;
plot(t_nl, theta1_nl, 'r-', 'LineWidth', 2);
yline(1, 'k:');
grid on;
legend('Linear Model', 'Non-Linear Model', 'Reference');
title('Part 7: Linear vs Non-Linear Step Response Comparison');
xlabel('Time (s)'); ylabel('Position (rad)');

% Output to Command Window
fprintf('--- PART 7: NON-LINEAR VALIDATION ---\n');
fprintf('Settling Time: %.4f s\n', Ts_nl);
fprintf('Overshoot:     %.2f %%\n\n', OS_nl);

%% --- PART 8: DISTURBANCE REJECTION (BOTH MODELS) ---
t_dist = 0:0.001:3;
x0_dist = [0; 0; 0; 0];
opts = odeset('RelTol', 1e-6, 'AbsTol', 1e-6);

% 1. Linear Model with Disturbance
[t_ld, x_ld] = ode45(@(t,x) linear_dist_dynamics(t, x, Kp, Kd, K_eq, J1, J2), t_dist, x0_dist, opts);
th1_ld = x_ld(:,1);

% 2. Non-Linear Model with Disturbance
[t_nd, x_nd] = ode45(@(t,x) nonlinear_dist_dynamics(t, x, Kp, Kd, k_lin, beta, J1, J2), t_dist, x0_dist, opts);
th1_nd = x_nd(:,1);

% Calculate Metrics
max_dev_lin = max(abs(th1_ld));
max_dev_nl = max(abs(th1_nd));

% Recovery Time (Threshold 0.001)
thresh = 0.001;
idx_rec_lin = find(abs(th1_ld) > thresh, 1, 'last');
if isempty(idx_rec_lin), rt_lin = 0; else, rt_lin = t_ld(idx_rec_lin); end

idx_rec_nl = find(abs(th1_nd) > thresh, 1, 'last');
if isempty(idx_rec_nl), rt_nl = 0; else, rt_nl = t_nd(idx_rec_nl); end

% FIGURE 3: PART 8
figure('Name', 'Part 8: Disturbance Rejection', 'Color', 'w');
subplot(2,1,1);
plot(t_ld, th1_ld, 'b--', 'LineWidth', 1.5); hold on;
plot(t_nd, th1_nd, 'r-', 'LineWidth', 1.5);
yline(thresh, 'k:', 'LineWidth', 1); yline(-thresh, 'k:', 'LineWidth', 1);
legend('Linear Model', 'Non-Linear Model');
title('Part 8: Motor Position under Arm Disturbance');
ylabel('Position (rad)'); grid on;

subplot(2,1,2);
dist_profile = arrayfun(@(t) get_disturbance(t), t_dist);
plot(t_dist, dist_profile, 'k', 'LineWidth', 1.5);
title('Applied Disturbance Torque');
xlabel('Time (s)'); ylabel('Torque (Nm)'); grid on; ylim([-0.5 2.5]);

% Output to Command Window
fprintf('--- PART 8: DISTURBANCE RESULTS ---\n');
fprintf('Linear Max Deviation: %.4f rad, Recovery: %.4f s\n', max_dev_lin, rt_lin);
fprintf('Non-Linear Max Dev:   %.4f rad, Recovery: %.4f s\n', max_dev_nl, rt_nl);


%% --- LOCAL FUNCTIONS ---

% A. Non-Linear Step Dynamics
function dxdt = nonlinear_dynamics(t, x, Kp, Kd, k, beta, J1, J2)
    th1 = x(1); dth1 = x(2); th2 = x(3); dth2 = x(4);
    ref = 1; 
    tau = Kp*(ref - th1) + Kd*(0 - dth1);
    dth = th1 - th2;
    T_spr = k*dth + beta*(dth)^3;
    dxdt = [dth1; (tau - T_spr)/J1; dth2; T_spr/J2];
end

% B. Linear Disturbance Dynamics
function dxdt = linear_dist_dynamics(t, x, Kp, Kd, K_eq, J1, J2)
    th1 = x(1); dth1 = x(2); th2 = x(3); dth2 = x(4);
    tau_m = -Kp*th1 - Kd*dth1;
    tau_d = get_disturbance(t);
    T_spr = K_eq * (th1 - th2);
    dxdt = [dth1; (tau_m - T_spr)/J1; dth2; (T_spr + tau_d)/J2];
end

% C. Non-Linear Disturbance Dynamics
function dxdt = nonlinear_dist_dynamics(t, x, Kp, Kd, k, beta, J1, J2)
    th1 = x(1); dth1 = x(2); th2 = x(3); dth2 = x(4);
    tau_m = -Kp*th1 - Kd*dth1;
    tau_d = get_disturbance(t);
    dth = th1 - th2;
    T_spr = k*dth + beta*(dth)^3;
    dxdt = [dth1; (tau_m - T_spr)/J1; dth2; (T_spr + tau_d)/J2];
end

% D. Disturbance Profile
function val = get_disturbance(t)
    if t<=0.25, val=8*t; elseif t<=0.5, val=2-8*(t-0.25); else, val=0; end
end