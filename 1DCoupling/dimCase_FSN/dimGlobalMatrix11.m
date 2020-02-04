function output = dimGlobalMatrix11(PoroProperty,input11)
%%             Comments 
% This function calculates the global matrix for solving the 1D HM
% problem in the following matrix form:
%
%   /                       \  /   \   /                           \ /    \  /                                  \
%   |[Mp]+ alpha*dt*[Kp],[C]|  |Pn |   |[Mp]+ (1-alpha)*dt*[Kp],[C]| |Pn_1|  |( alpha*{qn}+ (1-alpha)*{qn_1} )*dt|
%   |                       |* |   | = |                           |*|    |+ |                                   |
%   |-[C]^T            ,[Ku]|  |Un |   |          0           , 0  | |Un_1|  |           {Tn}                    |     
%   \                       /  \   /   \                           / \    /  \                                   /    


%==========================================================================
% 
%               Dimensional cases
% 
%==========================================================================
% for coupled hydraulic  problem:
%   
%   CM*dp/dt + b*dev/dt = k/mu * d2p/dx2
%
%   coeff_K = k/mu, coeff_M = CM, coeff_C = b;
%
% The discretized problem could be written as:
%
%  [Mp]*{P_dot} +[Cp]*{U_dot} + [Kp]*{P} + {q^d} = 0
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
%  -[Cp]^T*{P} + [Ku]*{U} + {T^d} = 0
%
%--------------------------------------------------------------------------

%% Input Parameters


k   = PoroProperty.k;   % intrinsic permeability [m2]
mu  = PoroProperty.mu;  % viscosity [Pa.s]
b   = PoroProperty.b;
CM  = PoroProperty.CM;  % storage coefficient [Pa^-1]
K   = PoroProperty.K;   % bulk modulus
G   = PoroProperty.G;   % shear modulus
lambda = K-2/3*G;

%% discretization
num_nodes   = input11.num_nodes;
L           = input11.L;
num_elem    = num_nodes -1;
x           = 0:L/(num_nodes-1):L;

input11.coord        = x;
input11.connec(:,1)  = 1:(num_nodes-1);
input11.connec(:,2)  = 2: num_nodes;


%% Global matrix for hydraulic problem


input11.mater(:,1) = k/mu *  ones(num_elem,1);
input11.mater(:,2) = CM   *  ones(num_elem,1);
input11.mater(:,3) = b    *  ones(num_elem,1);

[Kp,Mp,Cp]  = assembling11(input11);


%% Global matrix for poromechnanical problem
input11.mater(:,1) = ( lambda + 2*G ) *  ones(num_elem,1);
input11.mater(:,2) =        0         *  ones(num_elem,1);
input11.mater(:,3) =        b         *  ones(num_elem,1);

[Ku,~]  = assembling11(input11);

output.Kp = Kp;
output.Ku = Ku;
output.Mp = Mp;
output.Cp = Cp;
end 
