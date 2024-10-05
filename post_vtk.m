function [] = post_vtk(nodes,elements,type,q);
clc;
%% Renumbering nodes and elements
len_nodes = 1:length(nodes(:,1));
len_elements = 1:length(elements(:,1));
elements(:,1) = [];
for ii = 1:length(nodes(:,1));
B = find(elements(:,:)==nodes(ii,1));
elements(B) = len_nodes(ii);
endfor;
elements = [len_elements', elements];
nodes(:,1) = [];
nodes = [len_nodes', nodes];

%% Header information
filename = ['post_',type,'.vtk'];
version = '# vtk DataFile Version 3.0';
header = 'FEA Results';
format = 'ASCII';
geo_header = 'DATASET UNSTRUCTURED_GRID';

%% Geometry
if size(nodes,2) == 3;       % Adding a z-coordinates of 0 if 2D
nodes(:,4)=0;
else;
end;

switch(type);
case {"CPS9","CPS9R","CPE9","CPE9R"}
elements(:,end) = [];
end;

nodes(:,1) = [];             % Removing node numbers
num_nodes = numel(nodes)/3;  % Divide by 3, each cell has three coordinates 
points = ['POINTS ',num2str(num_nodes),' float'];

elements = elements(:,2:end);
num_ele = size(elements,1);
len_ele = size(elements,2);
cell_num(1:num_ele,1) = len_ele; % Length of cells
cells = [cell_num (elements-1)];
cell_size = numel(cells);
cell = ['CELLS ', num2str(num_ele),' ',num2str(cell_size)];

cell_type = ['CELL_TYPES ', num2str(num_ele)];
switch(type)
case {"CPS3","CPE3"}
cell_type2(1:num_ele,:) = 5;
case {"CPS4","CPE4","CPS4R"}
cell_type2(1:num_ele,:) = 9;      % Maybe 8? A pixel?
case {"CPS6","CPE6"}
cell_type2(1:num_ele,:) = 22;
case {"CPS8","CPE8"}
cell_type2(1:num_ele,:) = 23;
case {"CPS9","CPE9","CPS9R"}
cell_type2(1:num_ele,:) = 23;     % Middle node was removed above 
end


%% Post-processing results
points_dis = ['POINT_DATA ', num2str(num_nodes)];
scalar_type = 'SCALARS displacement double';
look_up = 'LOOKUP_TABLE default';       % Default look up table IE color gradient

q_x = q(1:2:end);
q_y = q(2:2:end);
q_z = zeros(size(q_x,1),1);
q_mag = sqrt(q_x.^2+q_y.^2);            % Magnitude of q at each node

deform_title = "VECTORS displacement double";
q_vect = [q_x q_y q_z];                 % Deformation of q at each node

%% Heading information
% Note the CD. Make a results folder
cd results
dlmwrite(filename,version,'delimiter','');
dlmwrite(filename,header,'-append','delimiter','');
dlmwrite(filename,format,'-append','delimiter','');

%% Geometry
dlmwrite(filename,geo_header,'-append','delimiter','');
dlmwrite(filename,points,'-append','delimiter','');
dlmwrite(filename,nodes,'-append','delimiter',' ','precision',3);
dlmwrite(filename,cell,'-append','delimiter','','roffset',1);
dlmwrite(filename,cells,'-append','delimiter',' ');
dlmwrite(filename,cell_type,'-append','delimiter','','roffset',1);
dlmwrite(filename,cell_type2,'-append');

%% Post-processing Results Gradient
dlmwrite(filename,points_dis,'-append','delimiter','','roffset',1);
dlmwrite(filename,scalar_type,'-append','delimiter','');
dlmwrite(filename,look_up,'-append','delimiter','');
dlmwrite(filename,q_mag,'-append','delimiter',' ','precision',4);

%%  Post-processing Results Deformation
dlmwrite(filename,deform_title,'-append','delimiter','','roffset',1);
dlmwrite(filename,q_vect,'-append','delimiter',' ','precision',4);
cd ..
endfunction;
