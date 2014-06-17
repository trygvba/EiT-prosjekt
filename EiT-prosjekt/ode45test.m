

ff =@(t,u)  [1,2,3;2,-4,3;1,0,4]*u;
u0=[0;1;2];
tspan=[0 5]

[tsol usol] =ode45(ff,tspan,u0);

eig([1,2,3;2,-4,3;1,0,4])

plot(tsol,usol)