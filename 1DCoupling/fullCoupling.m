function fullCoupling(input11)
% This function solves the coupled 1D HM problem using the sequential 
% method

%   /                       \  /   \   /                           \ /    \  /                                  \
%   |[Mp]+ alpha*dt*[Kp],[C]|  |Pn |   |[Mp]+ (1-alpha)*dt*[Kp],[C]| |Pn_1|  |( alpha*{qn}+ (1-alpha)*{qn_1} )*dt|
%   |                       |* |   | = |                           |*|    |+ |                                   |
%   |-[C]^T            ,[Ku]|  |Un |   |          0           , 0  | |Un_1|  |          -{Tn}                    |     
%   \                       /  \   /   \                           / \    /  \                                   /
% 
% which could be written as:
%
%   [mat1]*{vec_Xn} = [mat2]*{vec_Xn_1} + {vec_Yn} = {vec_Fn}
%                  
%-------------------------------------------------------------------------%
% input list:
% input11.
%         L           : total length of the domain
%         T           : total time for the study
%         num_iter    : number of iteration in one time step
%         num_tstep   : number of time step
%         num_nodes   : number of nodes 
%-------------------------------------------------------------------------%

alpha = 0.6;
L           = input11.L;
dt          = input11.T/input11.num_tstep;
num_nodes   = input11.num_nodes; 
num_iter    = input11.num_iter;

%% global matrix calculation
output      = globalMatrix11(PoroElasPara(),input11);
Mp          = output.Mp;
Kp          = output.Kp;
Cp          = output.Cp;
Ku          = output.Ku;

%% Loading parameters

q = zeros(num_nodes, num_tstep);
T = zeros(num_nodes, num_tstep);

%% initialize the unknowns
vec_X = zeros(2*num_nodes, num_tstep);


mat1                 = [Mp + alpha*dt*Kp, Cp ; - Cp', Ku ];
cond_mat1            = cond(mat1);
inv_mat1             = inv(mat1);
mat2                 = zeros(2*num_nodes,2*num_nodes);
mat2(1:num_nodes, :) = [Mp + (1-alpha)*dt*Kp, Cp];  
for n = 2 : num_tstep
    vec_Yn = [ (alpha*q(:,n) + (1-alpha)*q(:,n-1))*dt ; T(:,n)];
    vec_X(:,n) = inv_mat1 * ( mat2 *  vec_X(:,n-1) + vec_Yn );
end


end
