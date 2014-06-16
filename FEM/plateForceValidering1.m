function [ pf  ] = plateForceValidering1( f,loadrate,maxf,minf,t,period )


if t < period
    pf=0;
elseif (t >= period) && (t<2*period)
    pf=-period*loadrate + t*loadrate;
elseif (t>=2*period) && (t<3*period)
    pf=period*loadrate;
elseif (t>=3*period) && (t<4*period)
    pf=period*loadrate -(t-3*period)*loadrate;
else 
    pf=0;
end

end



