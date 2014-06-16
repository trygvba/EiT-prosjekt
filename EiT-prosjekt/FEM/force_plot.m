



maxf=10^-3;
min=200;

max=0;

f=zeros(1,1000);

t_max=10^-4;
dt=10^-8;


steps=t_max/dt;
p=t_max/9;

lr=maxf/p;

for i=1:steps
    t=i*dt;
    f(i)=plateForceValidering1(p,lr,max,min,t,p);
end

hold on
plot(f);