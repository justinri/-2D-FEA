function [kaa,kab,Force] = matrix_man(E,nu,nodes,elements,type,BC,force,idx,idx_BC);

%% Importing global matrix_man
%  Also making it a full matrix, no longer sparse
mat = assembly_select(E,nu,nodes,elements,type);

% Removes unneeded rows and columns
mat(idx,:) = [];
#mat(missing,:) = [];
kab = mat(:,idx_BC);

mat(:,idx) = [];
kaa = mat;


% Index force matrix
Force = zeros(2*max(nodes(:,1)),1);
BF = force(:,2) <= 1;
idx_force = [force(:,1)*2 - BF];
Force(idx_force) = force(:,3);
Force(idx) = [];
endfunction;

%{
  Old, slow code:
  That actually preforms matrix manipulation:
  (needs rows_below to pass through)
  % rows
rows_top = true(size(mat,1),1);
rows_top(rows_below) = false;
mod_mat = [mat(rows_top,:); mat(rows_below,:)];
len_row = 1:length(mat(rows_top,:)(:,1));

rows_below = [BC(:,1)*2 - B];


% columns
col_top = true(1,size(mod_mat,1));
col_top(rows_below) = false;
mod_mat = [mod_mat(:,col_top), mod_mat(:,rows_below)];
len_col = 1:length(mod_mat(:,col_top)(1,:));
len_col_ab = (length(mod_mat(:,col_top)(1,:))+1):length(mod_mat(1,:));

kaa = mod_mat(len_row,len_col);
kab = mod_mat(len_row,len_col_ab);
%}
