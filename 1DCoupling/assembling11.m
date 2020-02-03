function [K_glob,M_glob,C_glob]  = assembling11(input11)
%--------------------------------------------------------------------------
% function for assembling the global stiffness matrix for the 1D element

% input list:
% structure variable input11 with the following attributs:
% 
% coord: coordiate matrix
%
%   [x1 of glo_node_1; x1 of glo_node_2;...; x1 of glo_node_n]
%
%--------------------------------------------------------------------------
% connec: connectivity matrix relating the local node number to the global 
% node number
%
%           loc_node_i, loc_node_j
%   [glob_node_number_i; glob_node_number_j   -> element 1
%       ....                ....    
%    glob_node_number_i; glob_node_number_j;   -> element e
%       ....                ....
%    glob_node_number_i; glob_node_number_j]   -> element N_E
%
%--------------------------------------------------------------------------
% mater: material matrix
% 
%   [coeff_K, coeff_M; coeff_C;   -> element 1
%       ....          ....    
%    coeff_K, coeff_M; coeff_C;   -> element e
%       ....          ....
%    coeff_K, coeff_M; coeff_C]   -> element N_E
% 
%--------------------------------------------------------------------------
% for coupled hydraulic  problem:
%   
%   CM*dp/dt + b*dev/dt = k/mu * d2p/dx2
%
%   coeff_K = k/mu, coeff_M = CM, coeff_C = b;
%
% The discretized problem could be written as:
%
%  [M_glob]*{P_dot} +[C_glob]*{U_dot} + [K_glob_p]*{P} + {q^d} = 0
%
%--------------------------------------------------------------------------
% for coupled poromechanical problem:
%   
%   d( (lamda+2G)du/dx -b*p  )/dx =0
%   
%   coeff_K = lamda+2G, coeff_M = 0, coeff_C = b;
%
% The discretized problem could be written as:
%
%  -[C_glob]^T*{P} + [K_glob]*{U} + {T^d} = 0
%
%--------------------------------------------------------------------------

coord   = input11.coord;
connec  = input11.connect;
mater   = input11.mater;

num_Nodes = length(coord);
num_Elem  = length(connec);

K_glob = sparse(num_Nodes,num_Nodes);
M_glob = sparse(num_Nodes,num_Nodes);
C_glob = sparse(num_Nodes,num_Nodes);

elem_len = abs( connec(:,1)-connec(:,2) );

for e = 1:num_Elem
    [K_loc,M_loc,C_loc] = element11(elem_len(e));
    K_glob = K_glob(connec(e,:),connec(e,:)) + mater(e,1)* K_loc;
    M_glob = M_glob(connec(e,:),connec(e,:)) + mater(e,1)* M_loc;
    C_glob = C_glob(connec(e,:),connec(e,:)) + mater(e,1)* C_loc;
end



