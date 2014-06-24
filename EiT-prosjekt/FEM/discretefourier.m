function [F f] = discretefourier( Uz,steps,dt,MarkerNode)

t=(0:steps-1)*dt;
L=steps*dt;
Fs=1/dt;

NFFT = 2^nextpow2(steps);

F=fft(Uz(MarkerNode,:),NFFT)/L;

f = pi*Fs*linspace(0,1,NFFT/2+1);

%plot(f,abs(F(1:NFFT/2 +1)));
F=abs(F(1:NFFT/2 +1));


% s=sort(absF,'descend');
% mf = 1:50;
% for i = 1:50
%     mf(i) = 2*pi*f(absF==s(i));
% end
% subplot(2,1,2)
% plot(mf)
 
end

