function AU = AnalyticSolver(xx, tt, omega, wvel, ispulse, maxz)
plateDisp = @(t) -0.05*sin(omega*t).*(t>0).*(ispulse*t<(2*pi/omega));
AU = zeros(size(xx));
for i = 0:1:(10*wvel*tt+10)
    AU = AU + plateDisp(tt + (xx-(1+4*i)*maxz)/wvel) - plateDisp(tt-(xx+maxz*(3+4*i))/wvel);
end

end