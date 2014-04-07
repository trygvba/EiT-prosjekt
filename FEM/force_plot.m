



max=1000;
min=200;
lr=1600;


f=zeros(1,1000);

t_max=10^-4;
dt=10^-8;

steps=t_max/dt;
p=t_max/5;


for i=1:steps
    t=i*dt;
    f(i)=plateForceValidering1(p,lr,max,min,t,p);
end

plot(f);