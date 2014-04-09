function propagationPlotter(xplot)

r = [1 0 0];
g = [0 1 0];

N = size(xplot,1);
c = linspace(0,1,N);

figure
hold on;
for i=1:N
    plot(xplot(i,:),'color',c*r+(1-c)*g);
end
end