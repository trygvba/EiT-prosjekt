% OPTIMIZATION PLAYGROUND
% COPY CODE IN HERE AND RUN WITH PROFILER
I = eye(size(K2)); 
F_const = -Fmat_up*uz_up-Fmat_low*uz_low;
for i = 1:10
    steg1 = Mmod_inv*(F_cur+F_last);
    steg2 = (K1*Utemp_last);
    steg3 = -0.25*dt^2*(steg2 - steg1);
    steg4 = Utemp_last+dt*vel;
    steg5 = K2*(steg4 + steg3);
    X = K2*(Utemp_last-0.25*dt^2*(K1*Utemp_last));
   % Y = K2*((eye(size(K1))-(0.25*dt^2)*K1)*Utemp_last);
   
  
   
   
    Z = 2*(K2*Utemp_last) - Utemp_last;
    fitabrod = K2*(Utemp_last+dt*vel-0.25*dt^2*((K1*Utemp_last)-(Mmod_inv*(F_cur+F_last))));
    altfita =  K2*(2*Utemp_last + dt*vel + 0.25*dt^2*(Mmod_inv*(F_cur+F_last)))-Utemp_last;
    altaltfita = (I+(0.25*dt^2)*K1)\(2*Utemp_last + dt*vel + 0.25*dt^2*(Mmod_inv*(F_cur+F_last)))-Utemp_last;
end
sum(abs(altaltfita - altfita))
%sum(abs(X-Y))
%sum(abs(Y-Z))
sum(abs(X-Z))