function [ pF ] = plateForce2( t,omega,OLT )

  

if -omega*OLT*cos(omega*t) < 0
    pF=-1000;
else
    pF=0;
end

    


end

