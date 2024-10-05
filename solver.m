function [q] = solver(E,nu,nodes,elements,type,force,BC);
% Removing missing nodes numbers
BB = 1:max(nodes(:,1));
missingvalues = BB(~ismember(BB,nodes(:,1)));
idx_miss = [2*missingvalues-1; 2*missingvalues](:);
miss = zeros(2*size(nodes,1),1);
miss(idx_miss) = idx_miss;
missing = unique(miss);
missing(1) = [];
B = BC(:,2) <= 1;
idx_BC = [BC(1:end,1)*2 - B];
idx = sort([idx_BC; missing]);

% Importing kaa and kba
[kaa, kab, Force] = matrix_man(E,nu,nodes,elements,type,BC,force,idx,idx_BC);

% Solving for qa
qb = BC(:,3);
qa = kaa\(Force-kab*qb);


% putting all q's together
q = zeros(2*size(nodes,1),1);
q(idx_BC) = qb;
idx_qa = [1:2*max(nodes(:,1))].';
idx_qa(idx) = [];

q(idx_qa) = qa;


#whos idx_qa rows_below
#max(rows_below)
#rows_below(1:10)
#rows_below(335:end)
endfunction;

#% Removing missing nodes numbers
#BB = 1:max(nodes(:,1));
#missingvalues = BB(~ismember(BB,nodes(:,1)));
#idx_miss = [2*missingvalues-1; 2*missingvalues](:);
#miss = zeros(2*size(nodes,1),1);
#miss(idx_miss) = idx_miss;
#missing = unique(miss);
#missing(1) = [];

% Index force matrix
#Force = zeros(2*size(nodes,1),1);
#BF = force(:,2) <= 1;
#idx_force = [force(:,1)*2 - BF];
#Force(idx_force) = force(:,3);

#% Force manipulation, rows_below are BC locations
#rows_below = [BC(:,1)*2 - B];
#rows_top = true(size(Force,1),1);
#rows_top(rows_below) = false;
#force = [Force(rows_top,:)];


