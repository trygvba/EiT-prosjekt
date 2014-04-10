function AU = AnalyticSolver(xx, tt, omega, wvel)
%wvel = wvel/0.98;
plateDisp = @(t) -0.05*sin(omega*t).*(t>0).*(t<(2*pi/omega));
Base = @(x,t,n) (-1)^n*plateDisp(t-2*n/wvel + (((-1)^n)*x-1)/wvel);% - plateDisp(t-(x+3)/wvel).*(wvel*t>=(x+3)) + plateDisp(t+(x-5)/wvel).*(wvel*t>=(5-x)) - plateDisp(t-(x+7)/wvel).*(wvel*t>=(7+x));

AU = 0;
for i = 0:1:(floor(wvel*tt/2))
    AU = AU + Base(xx, tt, i);
end

end