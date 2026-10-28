clear; clc; close all;

J1 = 0.2;
J2 = 0.5;
k = 20;
beta = 100;
alpha0 = 0.1;
Keq = k + 3 * beta * alpha0^2;

w_res = sqrt(Keq*(J1+J2)/(J1*J2));
w0 = w_res;
zz = 0.1; zp = 0.7;
num_N = [1, 2*zz*w0, w0^2];
den_N = [1, 2*zp*w0, w0^2];
Notch = tf(num_N, den_N);

K = 3.5;
wz = 0.8; wp = 20;
Lead = K * tf([1/wz, 1], [1/wp, 1]);

C = Lead * Notch;

num_G1 = [J2, 0, Keq];
den_G1 = [J1*J2, 0, Keq*(J1+J2), 0, 0];
G1 = tf(num_G1, den_G1);
G2 = G1 * tf([Keq], [J2, 0, Keq]);

CL_th1 = feedback(C*G1, 1);
CL_th2 = minreal((G2*C) / (1 + G1*C));
CL_tau = minreal(C / (1 + C*G1));
CL_dist = minreal(G1 / (1 + C*G1)); 

t_sim = 0:0.001:5;

[y_lin_step, t_lin] = step(CL_th1, t_sim);

[num_C, den_C] = tfdata(C, 'v');
[Ac, Bc, Cc, Dc] = tf2ss(num_C, den_C);

x0 = zeros(4 + length(Ac), 1);
options = odeset('RelTol', 1e-6);

[t_nl_step, x_nl_step] = ode45(@(t,x) sys_dyn(t, x, J1, J2, k, beta, Ac, Bc, Cc, Dc, 1), t_sim, x0, options);
th1_nl_step = x_nl_step(:,1);
th2_nl_step = x_nl_step(:,3);
u_nl_step = zeros(length(t_nl_step), 1);
for i = 1:length(t_nl_step)
    xc = x_nl_step(i, 5:end)';
    err = 1 - th1_nl_step(i);
    u_nl_step(i) = Cc*xc + Dc*err;
end

[y_lin_dist, ~] = lsim(CL_dist, get_dist_profile(t_sim), t_sim);

[t_nl_dist, x_nl_dist] = ode45(@(t,x) sys_dyn(t, x, J1, J2, k, beta, Ac, Bc, Cc, Dc, 2), t_sim, x0, options);
th1_nl_dist = x_nl_dist(:,1);

info_lin_step = stepinfo(y_lin_step, t_lin, 1);
peak_lin_step = info_lin_step.Peak;
settle_lin_step = info_lin_step.SettlingTime;

peak_nl_step = max(th1_nl_step);
idx_s = find(abs(th1_nl_step - 1) > 0.02, 1, 'last');
if isempty(idx_s)
    settle_nl_step = 0;
else
    settle_nl_step = t_nl_step(idx_s);
end

peak_lin_dist = max(abs(y_lin_dist));
idx_ld = find(abs(y_lin_dist) > 0.02*peak_lin_dist, 1, 'last');
if isempty(idx_ld)
    rec_lin_dist = 0;
else
    rec_lin_dist = t_lin(idx_ld);
end

peak_nl_dist = max(abs(th1_nl_dist));
idx_nd = find(abs(th1_nl_dist) > 0.02*peak_nl_dist, 1, 'last');
if isempty(idx_nd)
    rec_nl_dist = 0;
else
    rec_nl_dist = t_nl_dist(idx_nd);
end

figure('Name', 'Part 4: Step Response with Metrics', 'Position', [100, 100, 800, 600]);
subplot(2,1,1);
plot(t_lin, y_lin_step, 'r--', 'LineWidth', 1.5); hold on;
plot(t_nl_step, th1_nl_step, 'b', 'LineWidth', 1.5);
yline(1, 'k:');
yline(1.02, 'g:'); yline(0.98, 'g:');
ylabel('\theta_1 (rad)'); title('Step Response: Motor Position'); 
legend('Linear', 'Nonlinear'); grid on;

str_step = { ...
    '\bfStep Response Metrics', ...
    '--------------------------------', ...
    sprintf('\\bfLinear Model:'), ...
    sprintf(' Peak: %.3f rad', peak_lin_step), ...
    sprintf(' Settling Time: %.3f s', settle_lin_step), ...
    ' ', ...
    sprintf('\\bfNonlinear Model:'), ...
    sprintf(' Peak: %.3f rad', peak_nl_step), ...
    sprintf(' Settling Time: %.3f s', settle_nl_step)};

dim_step = [0.55 0.5 0.3 0.3];
annotation('textbox', dim_step, 'String', str_step, 'FitBoxToText', 'on', 'BackgroundColor', 'white');

subplot(2,1,2);
plot(t_nl_step, u_nl_step, 'k', 'LineWidth', 1.2);
ylabel('Torque (Nm)'); title('Control Effort'); xlabel('Time (s)'); grid on;

figure('Name', 'Part 5: Disturbance Rejection with Metrics', 'Position', [150, 150, 800, 600]);
subplot(2,1,1);
plot(t_lin, get_dist_profile(t_lin), 'k', 'LineWidth', 1.5);
ylabel('Torque (Nm)'); title('Applied Disturbance Profile'); grid on;

subplot(2,1,2);
plot(t_lin, y_lin_dist, 'r--', 'LineWidth', 1.5); hold on;
plot(t_nl_dist, th1_nl_dist, 'b', 'LineWidth', 1.5);
ylabel('\theta_1 Deviation (rad)'); title('Disturbance Rejection Response'); 
xlabel('Time (s)'); legend('Linear', 'Nonlinear'); grid on;

str_dist = { ...
    '\bfDisturbance Metrics', ...
    '--------------------------------', ...
    sprintf('\\bfLinear Model:'), ...
    sprintf(' Max Deviation: %.4f rad', peak_lin_dist), ...
    sprintf(' Recovery Time: %.3f s', rec_lin_dist), ...
    ' ', ...
    sprintf('\\bfNonlinear Model:'), ...
    sprintf(' Max Deviation: %.4f rad', peak_nl_dist), ...
    sprintf(' Recovery Time: %.3f s', rec_nl_dist)};

dim_dist = [0.55 0.15 0.3 0.3];
annotation('textbox', dim_dist, 'String', str_dist, 'FitBoxToText', 'on', 'BackgroundColor', 'white');

function dxdt = sys_dyn(t, x, J1, J2, k, beta, Ac, Bc, Cc, Dc, mode)
    th1 = x(1); dth1 = x(2);
    th2 = x(3); dth2 = x(4);
    xc = x(5:end);
    
    if mode == 1
        ref = 1;
        td = 0;
    else
        ref = 0;
        if t < 0.25
            td = 8 * t;
        elseif t < 0.5
            td = 4 - 8 * t;
        else
            td = 0;
        end
    end
    
    err = ref - th1;
    tau = Cc*xc + Dc*err;
    
    dth = th1 - th2;
    Ts = k*dth + beta*(dth^3);
    
    ddth1 = (tau - Ts) / J1;
    ddth2 = (Ts + td) / J2;
    dxc = Ac*xc + Bc*err;
    
    dxdt = [dth1; ddth1; dth2; ddth2; dxc];
end

function dist = get_dist_profile(t)
    dist = zeros(size(t));
    for i = 1:length(t)
        tt = t(i);
        if tt < 0.25
            dist(i) = 8 * tt;
        elseif tt < 0.5
            dist(i) = 4 - 8 * tt;
        else
            dist(i) = 0;
        end
    end
end