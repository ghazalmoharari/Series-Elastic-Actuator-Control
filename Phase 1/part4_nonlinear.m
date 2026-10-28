function dydt = nonLinODE(t, y, J1, J2, k, beta, tau)
    theta1 = y(1);
    dtheta1 = y(2);
    theta2 = y(3);
    dtheta2 = y(4);
    
    if t >= 0
        T = tau;
    else
        T = 0;
    end
    
    delta = theta1 - theta2;
    T_spring = k * delta + beta * delta^3;
    
    ddtheta1 = (T - T_spring) / J1;
    ddtheta2 = T_spring / J2;
    
    dydt = [dtheta1; ddtheta1; dtheta2; ddtheta2];
end
J1 = 0.2;
J2 = 0.5;
k = 20;
beta = 100;
tau = 1;

y0 = [0; 0; 0; 0];
tspan = [0 5];

[t, y] = ode45(@(t, y) nonLinODE(t, y, J1, J2, k, beta, tau), tspan, y0);

theta1_nl = y(:, 1);
theta2_nl = y(:, 3);

figure;
plot(t, theta1_nl, 'b', 'LineWidth', 1.5);
hold on;
plot(t, theta2_nl, 'r', 'LineWidth', 1.5);
title('Non-Linear Model - Motor and Arm Angles');
xlabel('Time (s)');
ylabel('Angle (rad)');
legend('Motor Angle (θ1)', 'Arm Angle (θ2)');
grid on;
xlim([0 2]);


% Non-linear part (same as before)
J1 = 0.2;
J2 = 0.5;
k = 20;
beta = 100;
tau = 1;

y0 = [0; 0; 0; 0];
tspan = [0 5];

[t_nl, y_nl] = ode45(@(t, y) nonLinODE(t, y, J1, J2, k, beta, tau), tspan, y0);
theta1_nl = y_nl(:, 1);
theta2_nl = y_nl(:, 3);

% Linear part
fs = 20;
t_lin = 0 : 1/fs : 5 - 1/fs;
theta1_lin = 0.714*t_lin.^2 + 0.022*(1-cos(12.69*t_lin));
theta2_lin = 0.714*t_lin.^2 - 0.009*(1-cos(12.69*t_lin));

figure;

subplot(2,1,1);
plot(t_nl, theta1_nl, 'b', 'LineWidth', 1.5);
hold on;
plot(t_lin, theta1_lin, 'g--', 'LineWidth', 1.5);
title('Motor Angle (θ1) Comparison');
xlabel('Time (s)');
ylabel('Angle (rad)');
legend('Non-Linear', 'Linear');
grid on;
xlim([0 2]);

subplot(2,1,2);
plot(t_nl, theta2_nl, 'r', 'LineWidth', 1.5);
hold on;
plot(t_lin, theta2_lin, 'm--', 'LineWidth', 1.5);
title('Arm Angle (θ2) Comparison');
xlabel('Time (s)');
ylabel('Angle (rad)');
legend('Non-Linear', 'Linear');
grid on;
xlim([0 2]);