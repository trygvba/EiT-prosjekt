function Cp = StressMatrix(Ep, vp)
B1p = Ep/(2*vp.^2+vp-1)*[(vp-1) -vp -vp; -vp (vp-1) -vp; -vp -vp (vp-1)];
B2p = Ep/(2+2*vp)*eye(3);
Cp = [B1p zeros(3); zeros(3) B2p];
end