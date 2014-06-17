%Test script for the Newmark Beta-method:
%Problem is d^2y/dt^2=F(y,t);

%Let's try y"=-k*y
clear all;
close all;
k=0.01;

f = @(y) -k*y;

Beta = 0.25;

u0 = 1;
v0 = 0;

dt=0.1;
T = 200;
t=0:dt:T;
Nt=floor(T/dt);


v(1) = v0;
u(1) = u0;
acc(1) = f(u(1));
for i=2:Nt
    u(i) = (1/(1+Beta*dt^2*k))*(u(i-1)+dt*v(i-1)+((1-2*Beta)/2)*dt^2*acc(i-1));
    acc(i) = f(u(i));
    v(i) = v(i-1)+dt/2*(acc(i-1)+acc(i));
end
plot(t(1:(end-1)),u);