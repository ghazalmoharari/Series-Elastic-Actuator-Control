clear; clc; close all;

J1 = 0.2;
J2 = 0.5;
k_linear = 20;
beta = 100;
alpha_0 = 0.1;

K_eq = k_linear + 3 * beta * (alpha_0^2);

num_G = [J2, 0, K_eq];
den_G = [J1*J2, 0, K_eq*(J1+J2), 0, 0];
G = tf(num_G, den_G);

w_res = sqrt(K_eq * (J1 + J2) / (J1 * J2));

w0 = w_res;
zeta_z = 0.1;
zeta_p = 0.7;

num_Notch = [1, 2*zeta_z*w0, w0^2];
den_Notch = [1, 2*zeta_p*w0, w0^2];
N = tf(num_Notch, den_Notch);

w_z_lead = 0.8;
w_p_lead = 20.0;
K_gain = 3.5;

C_lead = K_gain * tf([1/w_z_lead, 1], [1/w_p_lead, 1]);

OL_Uncompensated = G;
OL_NotchOnly = N * G;
OL_Final = C_lead * N * G;

figure('Name', 'Bode Diagram Comparison');
opts = bodeoptions;
opts.FreqUnits = 'rad/s';
opts.PhaseMatching = 'on';
opts.Grid = 'on';
bode(OL_Uncompensated, opts);
hold on;
bode(OL_NotchOnly, opts);
bode(OL_Final, opts);
legend('Uncompensated (G)', 'With Notch (N*G)', 'Final (Lead*N*G)');
grid on;

[Gm, Pm, Wcg, Wcp] = margin(OL_Final);

fprintf('Gain Margin (Gm): %.2f dB\n', 20*log10(Gm));
fprintf('Phase Margin (Pm): %.2f deg\n', Pm);
fprintf('Gain Crossover Freq: %.2f rad/s\n', Wcg);
fprintf('Phase Crossover Freq: %.2f rad/s\n', Wcp);

if Pm > 45
    disp('RESULT: SUCCESS');
else
    disp('RESULT: FAIL');
end

CL_System = feedback(OL_Final, 1);
figure('Name', 'Closed Loop Step Response');
step(CL_System);
grid on;
stepinfo(CL_System)