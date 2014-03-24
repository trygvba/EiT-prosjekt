function [ pF ] = plateForce( t,omega )

  -1000*omega^2*sin(omega*t);

if -1000*omega^2*sin(omega*t) < 0
    pF=-1000*omega^2*sin(omega*t);
else
    pF=0;
end

    


end

