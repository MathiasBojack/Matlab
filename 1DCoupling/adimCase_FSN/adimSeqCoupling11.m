function [P,U,err] = adimSeqCoupling11(input11)
% This function solves the coupled 1D HM problem using the sequential 
% method

%   /                       \  /   \   /                           \ /    \  /                                  \
%   |[Mp]+ alpha*dt*[Kp],[C]|  |Pn |   |[Mp]+ (1-alpha)*dt*[Kp],[C]| |Pn_1|  |( alpha*{qn}+ (1-alpha)*{qn_1} )*dt|
%   |                       |* |   | = |                           |*|    |+ |                                   |
%   |-[C]^T            ,[Ku]|  |Un |   |          0           , 0  | |Un_1|  |           {Tn}                    |     
%   \                       /  \   /   \                           / \    /  \                                   /
% 
% which could be written as:
%
%   /           \ /   \ /          \ /    \  /         \   /         \
%   | K11 , K12 | |Pn | | L11, L12 | |Pn_1|  | {Fnp*dt}|   | {Y1n}   |
%   |           |*|   |=|          |*|    |+ |         | = |         |  
%   | K21 , K22 | |Un | |  0 , 0   | |Un_1|  |   {Tn}  |   | {Y2n}   |  
%   \           / \   / \          / \    /  \         /   \         /
%-------------------------------------------------------------------------%
% input list:
% input11.
%         L           : total length of the domain
%         T           : total time for the study
%         num_iter    : number of iteration in one time step
%         num_tstep   : number of time step
%         num_nodes   : number of nodes 
%-------------------------------------------------------------------------%
%%
alpha       = input11.alpha;
sigma0      = input11.sigma0;
dt          = input11.T/input11.num_tstep;
num_nodes   = input11.num_nodes; 
num_tstep   = input11.num_tstep;
num_iter    = input11.num_iter;
%% global matrix calculation
output      = adimGlobalMatrix11(PoroElasPara(),input11);
Mp          = output.Mp;
Kp          = output.Kp;
Cp          = output.Cp;
Ku          = output.Ku;

K11 =  Mp + alpha*dt*Kp;
K12 =  Cp;
K21 = -Cp';
K22 =  Ku;

L11 = Mp + (1-alpha)*dt*Kp;
L12 = Cp;

%% initialize the unknowns
P = zeros(num_nodes, num_tstep);
U = zeros(num_nodes, num_tstep);

%% Boundary conditions for mechanical problem


% fixed displacement at the bottom (X=0)
K21(1,:) = 0;
K22(1,:) = 0;
K22(1,1) = 1;
U(1,:)   = 0;


% drained at the top (X=1), pressure fixed at 0

K11(num_nodes, :) = 0;
K12(num_nodes, :) = 0;
K11(num_nodes,num_nodes) = 1;

L11(num_nodes, :) = 0;
L12(num_nodes, :) = 0;
L11(num_nodes,num_nodes) = 1;

P(num_nodes, :)   = 0;



%% Loading parameters

q = zeros(num_nodes, num_tstep);
T = zeros(num_nodes, num_tstep);
% no flow at X=0
q(1,:) = 0;
% applied stress sigma0 at X=1 for time step 2
T(num_nodes,2:num_tstep) = sigma0;


 err.pressure_first  = norm(full(K11^-1*K12*K22^-1*K21), 2);   
 err.disp_first  = norm(full(K22^-1*K21*K11^-1*K12), 2); 

for n = 2 : num_tstep
      Y1n =  L11 * P(:,n-1) + L12 * U(:,n-1) + dt*(alpha*q(:,n)+(1-alpha)*q(:,n-1));
      Y2n =  T(:,n);
      for i = 1:num_iter
            
            % ----------------------------------------------------------%
            % 1. start with the mechanical problem knowing the pressure %
            %-----------------------------------------------------------% 
            U(:,n) = K22\( Y2n - K21*P(:,n) ) ;
            P(:,n) = K11\( Y1n - K12*U(:,n) );

            % ------------------------------------------------------------%
            % 2. start with the hydraulic problem knowing the displacement%
            %-------------------------------------------------------------%
%             P(:,n) = K11\( Y1n - K12*U(:,n) )
%             U(:,n) = K22\( Y2n - K21*P(:,n) ) 
            
      end
end

end
