max_disp = load('../Data/maxdisp.mat');
top_disp = load('../Data/topdisp.mat');

max_disp = max_disp.max_disp'
top_disp = top_disp.topdisp;

v = 0.3;
E = 2*10^6;
R = 2.0267;

F = 4/3*E*sqrt(R)/(1-v^2)*max_disp.^(1.5);
 
figure
plot(max_disp,F,'r');
hold on
grid on
plot(top_disp,F,'b');
xlabel('Displacement')
ylabel('Force')
legend('Analytic','Simulated')


