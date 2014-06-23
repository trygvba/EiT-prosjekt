function [ max_freq, next_freq] = discretefourier( Uz,steps,dt,MarkerNode)

t=(0:steps-1)*dt;
L=steps*dt;
Fs=1/dt;

NFFT = 2^nextpow2(steps*dt)

F=fft(Uz(MarkerNode,:),NFFT)/L;

f = Fs/2*linspace(0,1,NFFT/2+1);


plot(f,abs(F(1:NFFT/2 +1)));
absF=abs(F(1:NFFT/2 +1));

max_freq=2*pi*f(absF==max(absF))
i1=find(absF==max(absF));
next_freq=2*pi*f((i1-1)+find(absF(i1:end)==max(absF(i1:end))));

end

