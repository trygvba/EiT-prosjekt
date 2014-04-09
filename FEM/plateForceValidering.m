function [ pf  ] = plateForceValidering( f,loadrate,maxf,minf,t,period )



if t < period
    pf=t*loadrate;
elseif (t>= period)&&(t<2*period)
    pf=period*loadrate
else
    pf=period*loadrate - (t-2*period)*loadrate;
end





end

