function [ pF  ] = plateForceValidering( f,loadrate,maxf,minf,t,period )



if cos(-pi/2 + pi/(2*period)*t) >= 0 
    pF=loadrate*t;
else
    pF=period*loadrate -loadrate*t;
end

    
    
    
    





end

