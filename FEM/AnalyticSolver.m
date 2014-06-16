function AU = AnalyticSolver(xx, tt, omega, wvel, ispulse)
plateDisp = @(t) -0.05*sin(omega*t).*(t>0).*(ispulse*t<(2*pi/omega));
AU = zeros(size(xx));
for i = 0:1:(10*wvel*tt+100)
    AU = AU + plateDisp((xx-1-4*i)/wvel+tt)- plateDisp(-(xx+3+4*i)/wvel+ tt);
end

end