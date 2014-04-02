function [ pf ] = plateForce3(t,omega,OLT,f  )

if cos(omega*t)>0
    pf=f*omega^2*OLT*cos(omega*t);
else
    pf=0;
end


end

