clear; clc; close all;

% --- 1. System Parameters ---
J1 = 0.2;       % Motor Inertia
J2 = 0.5;       % Arm Inertia
k = 20;         % Linear Stiffness
beta = 100;     % Nonlinear Stiffness Coefficient
alpha0 = 0.1;   % Operating Point

% Linear Equivalent Stiffness
Keq = k + 3 * beta * alpha0^2;

% --- 2. Linear Model Definitions ---
% Plant G1(s) = theta1(s) / tau(s)
num_G1 = [J2, 0, Keq];
den_G1 = [J1*J2, 0, Keq*(J1+J2), 0, 0];
G1 = tf(num_G1, den_G1);

% Plant G2(s) = theta2(s) / tau(s)
% Derived from: J2*s^2*theta2 = Keq*(theta1 - theta2) -> theta2/theta1 relation
% Easier: theta2 = (Keq / (J2*s^2 + Keq)) * theta1
% So G2 = G1 * (Keq / (J2*s^2 + Keq))
num_G2 = [Keq];
den_G2 = den_G1; % Same characteristic eq
G2 = tf(num_G2, den_G2);

% --- 3. Controller Design (Notch + Lead) ---
w_res = sqrt(Keq*(J1+J2)/(J1*J2));
w0 = w_res;
zz = 0.1; zp = 0.7;
Notch = tf([1, 2*zz*w0, w0^2], [1, 2*zp*w0, w0^2]);

K = 3.5;
wz = 0.8; wp = 20;
Lead = K * tf([1/wz, 1], [1/wp, 1]);

C = Lead * Notch;

% --- 4. Linear Simulation ---
% Closed Loop for Theta 1: theta1 = (GC)/(1+GC) * r
CL_theta1 = feedback(C*G1, 1);

% Closed Loop for Theta 2: theta2 = (G2*C)/(1+G1*C) * r
CL_theta2 = minreal((G2*C) / (1 + G1*C));

t = 0:0.01:5;
[y_lin_th1, t_lin] = step(CL_theta1, t);
[y_lin_th2, ~]     = step(CL_theta2, t);

% Calculate Linear Torque: tau = C * (r - theta1)
% Need transfer function Tau/R = C / (1 + C*G1)
CL_tau = minreal(C / (1 + C*G1));
[y_lin_tau, ~] = step(CL_tau, t);

% --- 5. Nonlinear Simulation ---
[num_C, den_C] = tfdata(C, 'v');
[Ac, Bc, Cc, Dc] = tf2ss(num_C, den_C);

x0 = zeros(4 + length(Ac), 1); % [th1, dth1, th2, dth2, xc...]
options = odeset('RelTol', 1e-6);
[t_nl, x_nl] = ode45(@(t,x) nonlinear_dyn(t, x, J1, J2, k, beta, Ac, Bc, Cc, Dc), t, x0, options);

theta1_nl = x_nl(:,1);
theta2_nl = x_nl(:,3);

% Reconstruct Nonlinear Torque for plotting
u_nl = zeros(length(t_nl), 1);
for i = 1:length(t_nl)
    xc = x_nl(i, 5:end)';
    err = 1 - theta1_nl(i);
    u_nl(i) = Cc*xc + Dc*err;
end

% --- 6. Plotting Angle Responses ---
figure('Name', 'Phase 2 Part 4: Angle Responses');

subplot(2,1,1);
plot(t_lin, y_lin_th1, 'r--', 'LineWidth', 1.5); hold on;
plot(t_nl, theta1_nl, 'b', 'LineWidth', 1.5);
ylabel('Angle (rad)');
legend('Linear Model', 'Nonlinear Model');
title('Motor Angle Response (\theta_1)');
grid on;

subplot(2,1,2);
plot(t_lin, y_lin_th2, 'r--', 'LineWidth', 1.5); hold on;
plot(t_nl, theta2_nl, 'b', 'LineWidth', 1.5);
ylabel('Angle (rad)'); xlabel('Time (s)');
legend('Linear Model', 'Nonlinear Model');
title('Arm Angle Response (\theta_2)');
grid on;

% Optional: Plot Torque (Required by item 4.1 text, usually)
figure('Name', 'Phase 2 Part 4: Torque Response');
plot(t_lin, y_lin_tau, 'r--', 'LineWidth', 1.5); hold on;
plot(t_nl, u_nl, 'b', 'LineWidth', 1.5);
ylabel('Torque (Nm)'); xlabel('Time (s)');
legend('Linear Model', 'Nonlinear Model');
title('Control Input Torque (\tau)');
grid on;

% --- Dynamics Function ---
function dxdt = nonlinear_dyn(t, x, J1, J2, k, beta, Ac, Bc, Cc, Dc)
    ref = 1; % Step Input
    
    % Unpack States
    theta1 = x(1); dtheta1 = x(2);
    theta2 = x(3); dtheta2 = x(4);
    xc = x(5:end);
    
    % Controller Output
    error = ref - theta1;
    tau = Cc*xc + Dc*error;
    
    % Nonlinear Spring
    dtheta = theta1 - theta2;
    T_spring = k*dtheta + beta*(dtheta^3);
    
    % Plant Dynamics
    ddtheta1 = (tau - T_spring) / J1;
    ddtheta2 = T_spring / J2;
    
    % Controller State Dynamics
    dxc = Ac*xc + Bc*error;
    
    dxdt = [dtheta1; ddtheta1; dtheta2; ddtheta2; dxc];
end