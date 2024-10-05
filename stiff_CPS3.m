function [k_glo_all] = stiff_CPS3(nodes,elements,Cs,k_glo_pre,k_glo_all);
for ii = 1:length(elements(:,1));
k_glo = k_glo_pre;
#
% Creating the x and y coordinates in to a vector
a1 = find(elements(ii,2)== nodes(:,1));
a2 = find(elements(ii,3)== nodes(:,1));
a3 = find(elements(ii,4)== nodes(:,1));
#
x=nodes([a1 a2 a3],2); 
y=nodes([a1 a2 a3],3);
#
x13 = x(1)-x(3);
x23 = x(2)-x(3);
x32 = -x23;
x21 = x(2)-x(1);
#
y13 = y(1)-y(3);
y23 = y(2)-y(3);
y31 = -y13;
y12 = y(1)-y(2);
#
#	n1 = xc;    % Shape functions are for reference
#	n2 = eta;
#	n3 = 1 - xc - eta;
#	Ja = ...
#	[x13, y13;
#	 x23, y23];
#    
DetJ = x13*y23 - x23*y13;
B1 = (1/DetJ)*[y23 0 y31 0 y12 0;
               0 x32 0 x13 0 x21;
	       x32 y23 x13 y31 x21 y12];
#
K1 = transpose(B1)*Cs*B1*(DetJ/2); % Where DetJ/2 = area
#
idx = [2*elements(ii,2:end)-1; 2*elements(ii,2:end)](:)';
k_glo([idx],[idx]) = K1;
k_glo_all = k_glo + k_glo_all;
endfor;
endfunction;
%{
 Old code for x and y coordinates:
 for jj = 2:1:length(elements(1,:));
a = nodes(:,1) == elements(ii,jj);
x(jj-1,1) = nodes(a,2);
y(jj-1,1) = nodes(a,3);
endfor;
%}
