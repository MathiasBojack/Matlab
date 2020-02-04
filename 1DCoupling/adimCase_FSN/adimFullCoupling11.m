function [P,U] = adimFullCoupling11(input11)
% This function solves the coupled 1D HM problem using the full coupling 
% method

%   /                       \  /   \   /                           \ /    \  /                                   \
%   |[Mp]+ alpha*dt*[Kp],[C]|  |Pn |   |[Mp]+ (1-alpha)*dt*[Kp],[C]| |Pn_1|  |( alpha*{qn}+ (1-alpha)*{qn_1} )*dt|
%   |                       |* |   | = |                           |*|    |+ |                                   |
%   |-[C]^T            ,[Ku]|  |Un |   |          0           , 0  | |Un_1|  |          {Tn}                     |     
%   \                       /  \   /   \                           / \    /  \                                   /
% 
% which could be written as:
%
%   [mat1]*{vec_Xn} = [mat2]*{vec_Xn_1} + {vec_Yn} = {vec_Fn}
%                  
%-------------------------------------------------------------------------%
% input list:
% input11.
%         T           : total time for the study
%         num_tstep   : number of time step
%         num_nodes   : number of nodes 
%         sigma0      : applied stress at the top 
%         alpha       : time discretization parameter
%-------------------------------------------------------------------------%
%%
alpha       = input11.alpha;
sigma0      = input11.sigma0;
dt          = input11.T/input11.num_tstep;
num_nodes   = input11.num_nodes; 
num_tstep   = input11.num_tstep;
%% global matrix calculation
output      = adimGlobalMatrix11(PoroElasPara(),input11);
Mp          = output.Mp;
Kp          = output.Kp;
Cp          = output.Cp;
Ku          = output.Ku;


mat1                 = [Mp + alpha*dt*Kp, Cp ; - Cp', Ku ];
mat2                 = zeros(2*num_nodes,2*num_nodes);
mat2(1:num_nodes, :) = [Mp + (1-alpha)*dt*Kp, Cp];  
%% initialize the unknowns
vec_X = zeros(2*num_nodes, num_tstep);

%% Boundary conditions for mechanical problem


% fixed displacement at the bottom (X=0)
mat1(num_nodes+1,:)           = 0;
mat1(num_nodes+1,num_nodes+1) = 1;

vec_X(num_nodes+1, :) = 0;

% drained at the top (X=1), pressure fixed at 0
mat1(num_nodes,:)              = 0;
mat1(num_nodes,num_nodes)      = 1;

mat2(num_nodes,:)              = 0;
mat2(num_nodes,num_nodes)      = 1;

vec_X(num_nodes, :) = 0;


%% Loading parameters

q = zeros(num_nodes, num_tstep );
T = zeros(num_nodes, num_tstep );
% no flow at X=0
q(1,:) = 0;
% applied stress sigma0 at X=1 for time step 2
T(num_nodes,2: num_tstep) = sigma0;

cond_mat1            = cond(full(mat1));


for n = 2 : num_tstep
    vec_Yn = [ (alpha*q(:,n) + (1-alpha)*q(:,n-1))*dt ; T(:,n)];
    vec_X(:,n) =  mat1\ ( mat2 *  vec_X(:,n-1) + vec_Yn );
end

P = vec_X(1:num_nodes,:);
U = vec_X(num_nodes+1:2*num_nodes,:);

end
