clc;
clear -x q1;
close all;

%% Material Constants
% Still need to add another layer, use Steel.E....
E = 200;%*10^3;  % MPa
nu = .25; 

% Type of elemnet to be used
type = "CPS3";

[nodes, elements, BC] = msh();
force = [1 1 100;
         1 2 100];
nodes(:,end) = [];

## Example Problem
#nodes = [1 5 5;
#         2 45 0;
#         3 65 45;
#         4 0 55;
#         5 25 2.5;
#         6 55 22.5;
#         7 32.5 50;
#         8 2.5 30];
##         9 28.75 26.25];
#elements = [1 1 2 3 4 5 6 7 8];
##            2 2 3 4];

#force = [1 1 100;
#         1 2 100];
#BC = [2 2 0];


% For problem 1
#[nodes, elements, BC, force] = mesh1(type);

 
% For problem 2
# [nodes, elements, BC, force] = mesh2(type);

% % For problem 3
 % fig = "right";
 % [nodes, elements, BC, force] = mesh3(fig);
 

 
% %% Displacement vector
tic;
q = solver(E,nu,nodes,elements,type,force,BC);
toc;
#q1 = q;
#isequal(q1,q)

% %% Post-processing in VTK format
post_vtk(nodes,elements,type,q);
