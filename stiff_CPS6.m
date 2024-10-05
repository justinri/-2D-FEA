function [K1] = stiff_CPS6(ii,nodes,elements,w,gp,Cs,K1);

% Creating the x and y coordinates in to a vector
x = nodes(elements(ii,2:end),2);
y = nodes(elements(ii,2:end),3);

for k = 2:2:8;
	xc = gp(k-1);
	eta = gp(k);
#	psi = gp(jj);
#	n1 = xc*(2*xc-1);
#	n2 = eta*(2*eta-1); % Shape functions are for reference
#	n3 = psi*(2*psi-1);
#	n4 = 4*xc*eta;
#	n5 = 4*psi*eta;
#	n6 = 4*xc*psi;
	Ja = ...
[ 4*x(6) - 3*x(3) - x(1) + 4*eta*x(3) + 4*eta*x(4) - 4*eta*x(5) - 4*eta*x(6) + 4*x(1)*xc + 4*x(3)*xc - 8*x(6)*xc, ...
  4*y(6) - 3*y(3) - y(1) + 4*eta*y(3) + 4*eta*y(4) - 4*eta*y(5) - 4*eta*y(6) + 4*xc*y(1) + 4*xc*y(3) - 8*xc*y(6);
  4*x(5) - 3*x(3) - x(2) + 4*eta*x(2) + 4*eta*x(3) - 8*eta*x(5) + 4*x(3)*xc + 4*x(4)*xc - 4*x(5)*xc - 4*x(6)*xc, ... 
  4*y(5) - 3*y(3) - y(2) + 4*eta*y(2) + 4*eta*y(3) - 8*eta*y(5) + 4*xc*y(3) + 4*xc*y(4) - 4*xc*y(5) - 4*xc*y(6)];

	DetJ = Ja(1,1)*Ja(2,2)-Ja(1,2)*Ja(2,1);
	
	A =1/DetJ*[Ja(2,2) -Ja(1,2) 0 0;
     		   0 0 -Ja(2,1) Ja(1,1);
	          -Ja(2,1) Ja(1,1) Ja(2,2) -Ja(1,2)];

	G = [(4*xc-1) 0 0 0 (4*eta+4*xc-3) 0 (4*eta) 0 ...
             (-4*eta) 0 (-4*eta-8*xc+4) 0;
             0 0 (4*eta-1) 0 (4*eta+4*xc-3) 0 (4*xc) 0 ...
             (-8*eta-4*xc+4) 0 (-4*xc) 0;
             0 (4*xc-1) 0 0 0 (4*eta+4*xc-3) 0 (4*eta) 0 ...
             (-4*eta) 0 (-4*eta-8*xc+4);
             0 0 0 (4*eta-1) 0 (4*eta+4*xc-3) 0 (4*xc) 0 ...
             (-8*eta-4*xc+4) 0 (-4*xc)];

	B1 = A*G;
	K = w(k/2)*transpose(B1)*Cs*B1*DetJ/2; % Where DetJ/2 = area
	K1 = K1 + K;
	
endfor;
endfunction;