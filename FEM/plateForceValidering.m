function [ pf  ] = plateForceValidering( f,loadrate,maxf,minf,t,period )



if t < period
    pf=t*loadrate;
elseif (t>= period)&&(t<2*period)
    pf=period*loadrate- (t-period)*loadrate;
elseif (t>=2*period) && (t<3*period)
    pf=(t-2*period)*loadrate;
else
    pf=period*loadrate - (t-3*period)*loadrate;
end






end

