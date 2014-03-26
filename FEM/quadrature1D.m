 % Approximates a 1-dimensional integral using Gaussian quadrature.
%Inputs: start- and endpoints (a,b), number of points Nq and the integrand g (as a function pointer).

function I = quadrature1D(a,b,Nq,g)
	y = @(x) (b-a)/2*x+(b+a)/2;
	switch Nq
		case 1
			p = 2;
			ztemp = 0;
		case 2
			p = [1 1];
			ztemp = (1/sqrt(3))*[-1 1];
		case 3
			p = (1/9)*[5 8 5];
			ztemp = sqrt(3/5)*[-1 0 1];
		case 4
			temp1 = 0.5;
			temp2 = sqrt(30)/36;
			p = [temp1-temp2 temp1+temp2 temp1+temp2 temp1-temp2];
			temp3 =3/7;
			temp4 = 2*sqrt(6/5)/7;
			ztemp = [-sqrt(temp3+temp4) -sqrt(temp3-temp4) sqrt(temp3-temp4) sqrt(temp3+temp4)];
		otherwise
			disp("Her ble det noe feil, er Nq = 1,2,3,4?")
			break;
	end
	z = feval(y,ztemp);
	I = (b-a)/2*feval(g,z)*p';
endfunction
