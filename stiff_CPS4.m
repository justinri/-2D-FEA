function [k_glo_all] = stiff_CPS4(nodes,elements,w,gp,Cs,k_glo_pre,k_glo_all);
   
for ii = 1:length(elements(:,1));
K1 = 0;
k_glo = k_glo_pre;

% Creating the x and y coordinates in to a vector
a1 = find(elements(ii,2)== nodes(:,1));
a2 = find(elements(ii,3)== nodes(:,1));
a3 = find(elements(ii,4)== nodes(:,1));
a4 = find(elements(ii,5)== nodes(:,1));
#loc = [a1 a2 a3];
x=nodes([a1 a2 a3 a4],2);
y=nodes([a1 a2 a3 a4],3);

for k = 1:2;
	for kk = 1:2;
	xc = gp(k);
	eta = gp(kk);
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
## w(k)*w(kk) Noting W's remove because they are 1, to save time.
	K = transpose(B1)*Cs*B1*DetJ;
	K1 = K1 + K;
	endfor;
endfor;

idx = [2*elements(ii,2:end)-1; 2*elements(ii,2:end)](:)';
k_glo([idx],[idx]) = sparse(K1);
k_glo_all = k_glo + k_glo_all;
   endfor;
endfunction;
