function [K1] = stiff_CPS4R(ii,nodes,elements,w,gp,Cs,K1);

% Creating the x and y coordinates in to a vector
x = nodes(elements(ii,2:end),2);
y = nodes(elements(ii,2:end),3);

	xc = gp;
	eta = gp;
#	n1 = .25*(1-xc)*(1-eta);  % Shape functions are for reference
#	n2 = .25*(1+xc)*(1-eta);
#	n3 = .25*(1+xc)*(1+eta);
#	n4 = .25*(1-xc)*(1+eta);
	Ja = .25*...
	[-(1-eta)*x(1)+(1-eta)*x(2)+(1+eta)*x(3)-(1+eta)*x(4),... 
	 -(1-eta)*y(1)+(1-eta)*y(2)+(1+eta)*y(3)-(1+eta)*y(4);
	 -(1-xc)*x(1)-(1+xc)*x(2)+(1+xc)*x(3)+(1-xc)*x(4),... 
	 -(1-xc)*y(1)-(1+xc)*y(2)+(1+xc)*y(3)+(1-xc)*y(4)];

	DetJ = Ja(1,1)*Ja(2,2)-Ja(1,2)*Ja(2,1);
	
	A =1/DetJ*[Ja(2,2) -Ja(1,2) 0 0;
     		   0 0 -Ja(2,1) Ja(1,1);
	          -Ja(2,1) Ja(1,1) Ja(2,2) -Ja(1,2)];

	G = .25*[-(1-eta) 0 1-eta 0 1+eta 0 -(1+eta) 0;
        	 -(1-xc) 0 -(1+xc) 0 (1+xc) 0 1-xc 0;
        	 0 -(1-eta) 0 1-eta 0 1+eta 0 -(1+eta);
         	 0 -(1-xc) 0 -(1+xc) 0 (1+xc) 0 1-xc];
	B1 = A*G;
	K1 = w*w*transpose(B1)*Cs*B1*DetJ;
endfunction;