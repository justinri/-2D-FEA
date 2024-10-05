function [K1] = stiff_CPS9(ii,nodes,elements,w,gp,Cs,K1);

% Creating the x and y coordinates in to a vector
x = nodes(elements(ii,2:end),2);
y = nodes(elements(ii,2:end),3);

loc = find(nodes(:,1)==elements(ii,2:end)');
x = nodes(loc,2);
y = nodes(loc,3);



[tf, loc] = ismember(nodes(:,1),elements(ii,2:end)');
tf
loc

for k = 1:3;
	for kk = 1:3;
	xc = gp(k);
	eta = gp(kk);
#	n1 = .25*eta*xc*(eta-1)*(xc-1); % Shape functions are for reference
#	n2 = .25*eta*xc*(eta-1)*(xc+1);
#	n3 = .25*eta*xc*(eta+1)*(xc+1);
#	n4 = .25*eta*xc*(eta+1)*(xc-1);
#	n5 = -.5*eta*(eta-1)*(xc-1)*(xc+1);
#	n6 = -.5*xc*(eta-1)*(eta+1)*(xc+1);
#	n7 = -.5*eta*(eta+1)*(xc-1)*(xc+1);
#	n8 = -.5*xc*(eta-1)*(eta+1)*(xc-1);
#	n9 = (eta-1)*(eta+1)*(xc-1)*(xc + 1);
	Ja = ...
[x(9)*(eta-1)*(eta+1)*(xc-1)-(x(8)*(eta-1)*(eta+1)*(xc-1))/2-(x(6)*(eta-1)*(eta+1)*(xc+1))/2+x(9)*(eta-1)*...
(eta+1)*(xc+1)+(eta*x(1)*xc*(eta-1))/4+(eta*x(2)*xc*(eta-1))/4+(eta*x(3)*xc*(eta+1))/4+(eta*x(4)*xc*(eta+1))/4+...
(eta*x(1)*(eta-1)*(xc-1))/4+(eta*x(2)*(eta-1)*(xc+1))/4+(eta*x(3)*(eta+1)*(xc+1))/4+(eta*x(4)*(eta+1)*(xc-1))/4-...
(eta*x(5)*(eta-1)*(xc-1))/2-(eta*x(5)*(eta-1)*(xc+1))/2-(x(6)*xc*(eta-1)*(eta+1))/2-(eta*x(7)*(eta+1)*(xc-1))/2-...
(eta*x(7)*(eta+1)*(xc+1))/2-(x(8)*xc*(eta-1)*(eta+1))/2, ...
%
y(9)*(eta-1)*(eta+1)*(xc-1)-(y(8)*(eta-1)*(eta+1)*(xc-1))/2-(y(6)*(eta-1)*(eta+1)*(xc+1))/2+y(9)*(eta-1)*(eta+1)*(xc+1)+...
(eta*xc*y(1)*(eta-1))/4+(eta*xc*y(2)*(eta-1))/4+(eta*xc*y(3)*(eta+1))/4+(eta*xc*y(4)*(eta+1))/4+(eta*y(1)*(eta-1)*(xc-1))/4+...
(eta*y(2)*(eta-1)*(xc+1))/4+(eta*y(3)*(eta+1)*(xc+1))/4+(eta*y(4)*(eta+1)*(xc-1))/4-(eta*y(5)*(eta-1)*(xc-1))/2-(eta*y(5)*...
(eta-1)*(xc+1))/2-(xc*y(6)*(eta-1)*(eta+1))/2-(eta*y(7)*(eta+1)*(xc-1))/2-(eta*y(7)*(eta+1)*(xc+1))/2-(xc*y(8)*(eta-1)*(eta+1))/2;
%
x(9)*(eta-1)*(xc-1)*(xc+1)-(x(7)*(eta+1)*(xc-1)*(xc+1))/2-(x(5)*(eta-1)*(xc-1)*(xc+1))/2+x(9)*(eta+1)*(xc-1)*(xc+1)+(eta*x(1)*xc...
*(xc-1))/4+(eta*x(2)*xc*(xc+1))/4+(eta*x(3)*xc*(xc+1))/4+(eta*x(4)*xc*(xc-1))/4+(x(1)*xc*(eta-1)*(xc-1))/4+(x(2)*xc*(eta-1)...
*(xc+1))/4+(x(3)*xc*(eta+1)*(xc+1))/4+(x(4)*xc*(eta+1)*(xc-1))/4-(eta*x(5)*(xc-1)*(xc+1))/2-(x(6)*xc*(eta-1)*(xc+1))/2-(x(6)...
*xc*(eta+1)*(xc+1))/2-(eta*x(7)*(xc-1)*(xc+1))/2-(x(8)*xc*(eta-1)*(xc-1))/2-(x(8)*xc*(eta+1)*(xc-1))/2, ...
%
y(9)*(eta-1)*(xc-1)*(xc+1)-(y(7)*(eta+1)*(xc-1)*(xc+1))/2-(y(5)*(eta-1)*(xc-1)*(xc+1))/2+y(9)*(eta+1)*(xc-1)*(xc+1)+...
(eta*xc*y(1)*(xc-1))/4+(eta*xc*y(2)*(xc+1))/4+(eta*xc*y(3)*(xc+1))/4+(eta*xc*y(4)*(xc-1))/4+(xc*y(1)*(eta-1)*(xc-1))/4+...
(xc*y(2)*(eta-1)*(xc+1))/4+(xc*y(3)*(eta+1)*(xc+1))/4+(xc*y(4)*(eta+1)*(xc-1))/4-(eta*y(5)*(xc-1)*(xc+1))/2-(xc*y(6)*(eta-1)...
*(xc+1))/2-(xc*y(6)*(eta+1)*(xc+1))/2-(eta*y(7)*(xc-1)*(xc+1))/2-(xc*y(8)*(eta-1)*(xc-1))/2-(xc*y(8)*(eta+1)*(xc-1))/2];

	DetJ = Ja(1,1)*Ja(2,2)-Ja(1,2)*Ja(2,1);
	
	A =1/DetJ*[Ja(2,2) -Ja(1,2) 0 0;
     		   0 0 -Ja(2,1) Ja(1,1);
	          -Ja(2,1) Ja(1,1) Ja(2,2) -Ja(1,2)];

	G = ...
[(eta*(2*xc-1)*(eta-1))/4,0,(eta*(2*xc+1)*(eta-1))/4,0,(eta*(2*xc+1)*(eta+1))/4,0,(eta*(2*xc-1)*(eta+1))/4,0,-eta*xc*(eta-1),0,...
-((eta^2-1)*(2*xc+1))/2,0,-eta*xc*(eta+1),0,-((eta^2-1)*(2*xc-1))/2,0,2*xc*(eta^2-1),0;
%
(xc*(2*eta-1)*(xc-1))/4,0,(xc*(2*eta-1)*(xc+1))/4,0,(xc*(2*eta+1)*(xc+1))/4,0,(xc*(2*eta+1)*(xc-1))/4,0,-((2*eta-1)*(xc^2-1))/2,...
0,-eta*xc*(xc+1),0,-((2*eta+1)*(xc^2-1))/2,0,-eta*xc*(xc-1),0,2*eta*(xc^2-1),0;
%
0,(eta*(2*xc-1)*(eta-1))/4,0,(eta*(2*xc+1)*(eta-1))/4,0,(eta*(2*xc+1)*(eta+1))/4,0,(eta*(2*xc-1)*(eta+1))/4,0,-eta*xc*(eta-1),...
0,-((eta^2-1)*(2*xc+1))/2,0,-eta*xc*(eta+1),0,-((eta^2-1)*(2*xc-1))/2,0,2*xc*(eta^2-1);
%
0,(xc*(2*eta-1)*(xc-1))/4,0,(xc*(2*eta-1)*(xc+1))/4,0,(xc*(2*eta+1)*(xc+1))/4,0,(xc*(2*eta+1)*(xc-1))/4,0,-((2*eta-1)*(xc^2-1))/2,...
0,-eta*xc*(xc+1),0,-((2*eta+1)*(xc^2-1))/2,0,-eta*xc*(xc-1),0,2*eta*(xc^2-1)];

whos A G
adad
	B1 = A*G;
	K = w(k)*w(kk)*transpose(B1)*Cs*B1*DetJ;
	K1 = K1 + K;
	endfor;
endfor;

endfunction;

