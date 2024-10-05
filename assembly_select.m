function [k_glo_all] = assembly_select(E,nu,nodes,elements,type);
% Pre-allocating Space
m_size = 2*max(nodes(:,1));
size_ele = (2*length(elements(1,2:end)))^2;
rittenhouse_constant = .4;
size_node = size_ele*size(elements,1)*rittenhouse_constant;
k_glo_pre = spalloc(m_size,m_size,size_ele);
k_glo_all = spalloc(m_size,m_size,size_node);

%{ 
  Picks the correct element type.
  Creates stiffness matrix for each element
  Storing each in it's own cell
%}
switch(type)
%% Plane Stress
	case {"CPS3", "CPS6", "CPS4","CPS4R", "CPS8", "CPS9","CPS9R"}
	Cs = E/(1-nu^2)*[1 nu 0;
		         nu 1 0;
			 0 0 (1-nu)/2];
   switch(type)
   case "CPS3" % Plane Stress Triangle
   k_glo_all = stiff_CPS3(nodes,elements,Cs,k_glo_pre,k_glo_all);
   
   case "CPS6" % Plane Stress (Second Order) Triangle
   w = [-27/48 25/48 25/48 25/48];
   gp =[1/3 1/3 3/5 1/5 1/5 1/5 1/5 3/5];
   for ii = 1:length(elements(:,1));
   K1 = 0;
   k_glo = k_glo_pre;
   K1 = stiff_CPS6(ii,nodes,elements,w,gp,Cs,K1);
   idx = [2*elements(ii,2:end)-1; 2*elements(ii,2:end)](:)';
   k_glo([idx],[idx]) = sparse(K1);
   k_glo_all = k_glo + k_glo_all;
   endfor;
   
   case "CPS4" % Plane Stress Quad
   w = [1 1];
   gp = [-1/sqrt(3) 1/sqrt(3)];
   k_glo_all = stiff_CPS4(nodes,elements,w,gp,Cs,k_glo_pre,k_glo_all);
   
   case "CPS4R" % Plane Stress Quad Reduced Integration
   w = 2;
   gp = 0;
   for ii = 1:length(elements(:,1));
   K1 = 0;
   k_glo = k_glo_pre;
   K1 = stiff_CPS4R(ii,nodes,elements,w,gp,Cs,K1);
   idx = [2*elements(ii,2:end)-1; 2*elements(ii,2:end)](:)';
   k_glo([idx],[idx]) = sparse(K1);
   k_glo_all = k_glo + k_glo_all;
   endfor;
   
   case "CPS8" % Plane Stress Quad-Serendipity
   w = [5/9 5/9 8/9];
   gp = [sqrt(3/5) -sqrt(3/5) 0];
   k_glo_all = stiff_CPS8(nodes,elements,w,gp,Cs,k_glo_pre,k_glo_all);
   
   case "CPS9" % Plane Stress (Second Order) Quad
   w = [5/9 5/9 8/9];
   gp = [sqrt(3/5) -sqrt(3/5) 0];
   for ii = 1:length(elements(:,1));
   K1 = 0;
   k_glo = k_glo_pre;
   K1 = stiff_CPS9(ii,nodes,elements,w,gp,Cs,K1);
   idx = [2*elements(ii,2:end)-1; 2*elements(ii,2:end)](:)';
   k_glo([idx],[idx]) = sparse(K1);
   k_glo_all = k_glo_all + k_glo;
   endfor;
   
   case "CPS9R" % Plane Stress (Second Order) Quad Reduced Integration
   w = [1 1];
   gp = [-1/sqrt(3) 1/sqrt(3)];
   for ii = 1:length(elements(:,1));
   K1 = 0;
   k_glo = k_glo_pre;
   K1 = stiff_CPS9R(ii,nodes,elements,w,gp,Cs,K1);
   idx = [2*elements(ii,2:end)-1; 2*elements(ii,2:end)](:)';
   k_glo([idx],[idx]) = sparse(K1);
   k_glo_all = k_glo + k_glo_all;
   endfor;
   end;   
      
   
 %% Plane Strain
 	case {"CPE3", "CPE6", "CPE4","CPE4R", "CPE8", "CPE9","CPE9R"}
Cs = E/((1+nu)*(1-2*nu))*[(1-nu) nu 0;
                          nu (1-nu) 0;
                          0 0 (1-2*nu)];
 switch(type)
   case "CPE3" % Plane Strain Triangle
   for ii = 1:length(elements(:,1));
   k_glo = k_glo_pre;
   K1 = stiff_CPS3(ii,nodes,elements,Cs);
   idx = [2*elements(ii,2:end)-1; 2*elements(ii,2:end)](:)';
   k_glo([idx],[idx]) = K1;
   k_glo_all{ii} = k_glo;
   endfor;
      
   case "CPE6" % Plane Strain (Second Order) Triangle
   w = [-27/48 25/48 25/48 25/48];
   gp =[1/3 1/3 3/5 1/5 1/5 1/5 1/5 3/5];
   for ii = 1:length(elements(:,1));
   K1 = 0;
   k_glo = k_glo_pre;
   K1 = stiff_CPS6(ii,nodes,elements,w,gp,Cs,K1);
   idx = [2*elements(ii,2:end)-1; 2*elements(ii,2:end)](:)';
   k_glo([idx],[idx]) = K1;
   k_glo_all{ii} = k_glo;
   endfor;
   
   case "CPE4" % Plane Strain Quad
   w = [1 1];
   gp = [-1/sqrt(3) 1/sqrt(3)];
   for ii = 1:length(elements(:,1));
   K1 = 0;
   k_glo = k_glo_pre;
   K1 = stiff_CPS4(ii,nodes,elements,w,gp,Cs,K1);
   idx = [2*elements(ii,2:end)-1; 2*elements(ii,2:end)](:)';
   k_glo([idx],[idx]) = K1;
   k_glo_all{ii} = k_glo;
   endfor;
   
   case "CPE4R" % Plane Strain Quad Reduced Integration
   w = 2;
   gp = 0;
   for ii = 1:length(elements(:,1));
   K1 = 0;
   k_glo = k_glo_pre;
   K1 = stiff_CPS4R(ii,nodes,elements,w,gp,Cs,K1);
   idx = [2*elements(ii,2:end)-1; 2*elements(ii,2:end)](:)';
   k_glo([idx],[idx]) = K1;
   k_glo_all{ii} = k_glo;
   endfor;
   
   case "CPE8" % Plane Strain Quad-Serendipity
   w = [5/9 5/9 8/9];
   gp = [sqrt(3/5) -sqrt(3/5) 0];
   for ii = 1:length(elements(:,1));
   K1 = 0;
   k_glo = k_glo_pre;
   K1 = stiff_CPS8(ii,nodes,elements,w,gp,Cs,K1);
   idx = [2*elements(ii,2:end)-1; 2*elements(ii,2:end)](:)';
   k_glo([idx],[idx]) = K1;
   k_glo_all{ii} = k_glo;
   endfor;
   
   case "CPE9" % Plane Strain (Second Order) Quad
   w = [5/9 5/9 8/9];
   gp = [sqrt(3/5) -sqrt(3/5) 0];
   for ii = 1:length(elements(:,1));
   K1 = 0;
   k_glo = k_glo_pre;
   K1 = stiff_CPS9(ii,nodes,elements,w,gp,Cs,K1);
   idx = [2*elements(ii,2:end)-1; 2*elements(ii,2:end)](:)';
   k_glo([idx],[idx]) = K1;
   k_glo_all{ii} = k_glo;
   endfor;
   
   case "CPE9R" % Plane Strain (Second Order) Quad Reduced Integration
   w = [1 1];
   gp = [-1/sqrt(3) 1/sqrt(3)];
   for ii = 1:length(elements(:,1));
   K1 = 0;
   k_glo = k_glo_pre;
   K1 = stiff_CPS9R(ii,nodes,elements,w,gp,Cs,K1);
   idx = [2*elements(ii,2:end)-1; 2*elements(ii,2:end)](:)';
   k_glo([idx],[idx]) = K1;
   k_glo_all{ii} = k_glo;
   endfor;
 
end;
end;
endfunction;
