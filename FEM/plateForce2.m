function [ pF ] = plateForce2( t,omega,OLT,f )

  

if -omega*OLT*cos(omega*t) < 0
    pF=-f;
else
    pF=0;
end

    


end

