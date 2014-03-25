function [ pF ] = plateForce( t,omega,OLT )

  

if -OLT*omega^2*sin(omega*t) < 0
    pF=-OLT*omega^2*sin(omega*t);
else
    pF=0;
end

    


end

