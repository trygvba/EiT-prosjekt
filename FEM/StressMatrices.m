function [Cp Cs] = StressMatrices(Ep, Es, vp, vs)
B1p = Ep/(2*vp^2+vp-1)*[(vp-1) -vp -vp; -vp (vp-1) -vp; -vp -vp (vp-1)];
B2p = Ep/(2+2*vp)*eye(3);
Cp = [B1p zeros(3); zeros(3) B2p];

B1s = Es/(2*vs^2+vs-1)*[(vs-1) -vs -vs; -vs (vs-1) -vs; -vs -vs (vs-1)];
B2s = Es/(2+2*vs)*eye(3);
Cs = [B1s zeros(3); zeros(3) B2s];

end

