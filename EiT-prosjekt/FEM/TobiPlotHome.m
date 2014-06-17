%TobiPlot
node2 = 779;
xx = p(node2,3);

norm = @(vector, LP) (sum(abs(vector).^LP)).^(1/LP);

Asol = zeros(size(Uz));
h = waitbar(0, 'Solving analytically');
for tue = 1:steps
    Asol(:,tue) = AnalyticSolver(p(:,3),tue*dt,omega,wvel, ispulse, 1);
    waitbar(tue/steps);
end
close(h);

Lp = 2;
error = max(abs(Uz - Asol)/0.05);
iterr = zeros(steps,1);
iterr1 = iterr; 
cumerr = iterr;
for i = 1:steps
    iterr(i) = norm(error(1:i),Lp)/i;
    cumerr(i) = norm(error(1:i), Lp)/steps;
%    iterr1(i) = norm(error(1:i),1)/i;
end

figHan = figure;
set(figHan, 'Position', [0, 0, 500, 500]);
subplot(2,1,1)
hold on
plot(1:steps,Uz(node2,:))
plot(1:steps,Asol(node2,:), 'r');

%plot(asol*wvel, Uz(node2,:)-Asol, 'g')
legend('Simulation', 'Analytic');
title(['z_0 = ' num2str(p(node2,3)) ', n = ' num2str(harmonicMode)]);
ylabel('Displacement');
subplot(2,1,2)
hold on

plot(iterr);
plot(cumerr,'r')
%plot(iterr1,'g');
title(['l^\infty(x,y,z) \times l^' num2str(Lp) '(t) norm error']);
legend('Iterative error','Cumulative error')%,'Absolute error'
xlabel('Timesteps');
ylabel('Error percentage');
disp(cumerr(end))