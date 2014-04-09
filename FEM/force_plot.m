



maxf=10^-3;
min=200;



f=zeros(1,1000);

t_max=10^-4;
dt=10^-8;


steps=t_max/dt;
p=t_max/3;

lr=maxf/p;

for i=1:steps
    t=i*dt;
    f(i)=plateForceValidering(p,lr,max,min,t,p);
end

hold on
plot(f);