function [k_glo_all] = stiff_CPS8(nodes,elements,w,gp,Cs,k_glo_pre,k_glo_all);
for ii = 1:length(elements(:,1));
K1 = 0;
k_glo = k_glo_pre;

% Creating the x and y coordinates in to a vector
a1 = find(elements(ii,2)== nodes(:,1));
a2 = find(elements(ii,3)== nodes(:,1));
a3 = find(elements(ii,4)== nodes(:,1));
a4 = find(elements(ii,5)== nodes(:,1));
a5 = find(elements(ii,6)== nodes(:,1));
a6 = find(elements(ii,7)== nodes(:,1));
a7 = find(elements(ii,8)== nodes(:,1));
a8 = find(elements(ii,9)== nodes(:,1));
#loc = [a1 a2 a3];
x=nodes([a1 a2 a3 a4 a5 a6 a7 a8],2);
y=nodes([a1 a2 a3 a4 a5 a6 a7 a8],3);

for k = 1:3;
	for kk = 1:3;
	xc = gp(k);
	eta = gp(kk);
#   n1 = -.25*(1-xc)*(1-eta)*(1+xc+eta); % Shape functions are for reference
#	n2 = -.25*(1+xc)*(1-eta)*(1-xc+eta);
#	n3 = -.25*(1+xc)*(1+eta)*(1-xc-eta);
#	n4 = -.25*(1-xc)*(1+eta)*(1+xc-eta);
#	n5 = .5*(1-xc^2)*(1-eta);
#	n6 = .5*(1+xc)*(1-eta^2);
#	n7 = .5*(1-xc^2)*(1+eta);
#	n8 = .5*(1-xc)*(1-eta^2);

	Ja = ...
[(x(8)*(eta^2-1))/2-(x(6)*(eta^2-1))/2+x(5)*xc*(eta-1)-x(7)*xc*(eta+1)+(x(2)*(eta-1)*(eta-xc+1))/4+(x(4)*(eta+1)* ...
(xc-eta+1))/4-x(1)*(xc/4-1/4)*(eta-1)-x(2)*(xc/4+1/4)*(eta-1)+x(3)*(xc/4+1/4)*(eta+1)+x(4)*(xc/4-1/4)*(eta+1)-...
(x(1)*(eta-1)*(eta+xc+1))/4+(x(3)*(eta+1)*(eta+xc-1))/4, ...
#
(y(8)*(eta^2-1))/2-(y(6)*(eta^2-1))/2+xc*y(5)*(eta-1)-xc*y(7)*(eta+1)+(y(2)*(eta-1)*(eta-xc+1))/4+(y(4)*(eta+1)...
*(xc-eta+1))/4-y(1)*(xc/4-1/4)*(eta-1)-y(2)*(xc/4+1/4)*(eta-1)+y(3)*(xc/4+1/4)*(eta+1)+y(4)*(xc/4-1/4)*(eta+1)-...
(y(1)*(eta-1)*(eta+xc+1))/4+(y(3)*(eta+1)*(eta+xc-1))/4;
#
x(5)*(xc^2/2-1/2)-x(7)*(xc^2/2-1/2)-x(1)*(xc/4-1/4)*(eta+xc+1)+x(3)*(xc/4+1/4)*(eta+xc-1)-2*eta*x(6)*(xc/2+1/2)+...
2*eta*x(8)*(xc/2-1/2)+x(2)*(xc/4+1/4)*(eta-xc+1)+x(4)*(xc/4-1/4)*(xc-eta+1)-x(1)*(xc/4-1/4)*(eta-1)+x(2)*(xc/4+1/4)...
*(eta-1)+x(3)*(xc/4+1/4)*(eta+1)-x(4)*(xc/4-1/4)*(eta+1), ...
#
y(5)*(xc^2/2-1/2)-y(7)*(xc^2/2-1/2)-y(1)*(xc/4-1/4)*(eta+xc+1)+y(3)*(xc/4+1/4)*(eta+xc-1)-2*eta*y(6)*(xc/2+1/2)+2*eta*...
y(8)*(xc/2-1/2)+y(2)*(xc/4+1/4)*(eta-xc+1)+y(4)*(xc/4-1/4)*(xc-eta+1)-y(1)*(xc/4-1/4)*(eta-1)+y(2)*(xc/4+1/4)*(eta-1)+...
y(3)*(xc/4+1/4)*(eta+1)-y(4)*(xc/4-1/4)*(eta+1)];

 

	DetJ = Ja(1,1)*Ja(2,2)-Ja(1,2)*Ja(2,1);
	
	A =1/DetJ*[Ja(2,2) -Ja(1,2) 0 0;
     		   0 0 -Ja(2,1) Ja(1,1);
	          -Ja(2,1) Ja(1,1) Ja(2,2) -Ja(1,2)];

		G=...
[-((eta+2*xc)*(eta-1))/4,0,((eta-2*xc)*(eta-1))/4,0,((eta+2*xc)*(eta+1))/4,0,-((eta-2*xc)*(eta+1))/4,0,xc*(eta-1),0,...
1/2-eta^2/2,0,-xc*(eta+1),0,eta^2/2-1/2,0;
#
-((2*eta+xc)*(xc-1))/4,0,((xc+1)*(2*eta-xc))/4,0,((2*eta+xc)*(xc+1))/4,0,-((xc-1)*(2*eta-xc))/4,0,xc^2/2-1/2,0,...
-eta*(xc+1),0,1/2-xc^2/2,0,eta*(xc-1),0;
#
0,-((eta+2*xc)*(eta-1))/4,0,((eta-2*xc)*(eta-1))/4,0,((eta+2*xc)*(eta+1))/4,0,-((eta-2*xc)*(eta+1))/4,0,xc*(eta-1),...
0,1/2-eta^2/2,0,-xc*(eta+1),0,eta^2/2-1/2;
#
0,-((2*eta+xc)*(xc-1))/4,0,((xc+1)*(2*eta-xc))/4,0,((2*eta+xc)*(xc+1))/4,0,-((xc-1)*(2*eta-xc))/4,0,xc^2/2-1/2,0,...
-eta*(xc+1),0,1/2-xc^2/2,0,eta*(xc-1)];


	B1 = A*G;
	K = w(k)*w(kk)*transpose(B1)*Cs*B1*DetJ;
	K1 = K1 + K;
	endfor;
endfor;
idx = [2*elements(ii,2:end)-1; 2*elements(ii,2:end)](:)';
k_glo([idx],[idx]) = sparse(K1);
k_glo_all = k_glo + k_glo_all;
endfor;
endfunction;
