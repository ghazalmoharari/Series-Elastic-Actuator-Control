J1=0.2;
J2=0.5;
k=20;
beta=100;
alpha0=0.1;
Keq=k*alpha0+beta*alpha0^3;
A=[0 1 0 0; -115 0 115 0;0 0 0 1;46 0 -46 0 ];
B=[0;5;0;0];
C=[1 0 0 0; 0 0 1 0];
D=[0;0];
[num, denum]=ss2tf(A,B,C,D);
T = 5;
fs = 20;
t = 0 : 1/fs : T - 1 / fs;

teta1 = 0.714*t.^2+0.022*(1-cos(12.69*t));
teta2=  0.714*t.^2-0.009*(1-cos(12.69*t));
teta1Plot=plot(t,teta1, LineWidth=2, Color='GREEN');
hold on
teta2Plot=plot(t,teta2, LineWidth=2, Color='RED', LineStyle='--');
title(' Linearized plots of teta1(t) and teta2(t) ');
xlabel("Time(s)");
ylabel("Values:")
legend([teta1Plot, teta2Plot], {'teta1(t)', 'teta2(t)'})

xlim([0 2]);
ylim([0 5]);

