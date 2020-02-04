function output = adimGlobalMatrix11(PoroProperty,input11)
%%             Comments 
% This function calculates the global matrix for solving the 1D HM
% problem in the following matrix form:
%
%   /                       \  /   \   /                           \ /    \  /                                  \
%   |[Mp]+ alpha*dt*[Kp],[C]|  |Pn |   |[Mp]+ (1-alpha)*dt*[Kp],[C]| |Pn_1|  |( alpha*{qn}+ (1-alpha)*{qn_1} )*dt|
%   |                       |* |   | = |                           |*|    |+ |                                   |
%   |-[C]^T            ,[Ku]|  |Un |   |          0           , 0  | |Un_1|  |          -{Tn}                    |     
%   \                       /  \   /   \                           / \    /  \                                   /    

%==========================================================================
% 
%               1. Adimensional cases
% 
%==========================================================================
% for coupled hydraulic  problem: Fixed strain split
%   
%   dP/dT + b*dev/dT = d2P/dX2
% 
%   where the adimensional parameters are defined by
%   
%   P = p*CM, X = x/L, T = t/T0, T0 = mu*CM*L^2/k
%
%   coeff_K = 1, coeff_M = 1, coeff_C = b;
%
% The discretized problem could be written as:
%
%  [Mp]*{P_dot} +[Cp]*{U_dot} + [Kp]*{P} + {q^d} = 0
%
%--------------------------------------------------------------------------
% for coupled poromechanical problem:
%   
%   d( (lamda+2G)*CM)*dU/dX -b*P  )/dX =0
%   
%   coeff_K = (lamda+2G)*CM, coeff_M = 0, coeff_C = b;
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

%% discretization
num_nodes   = input11.num_nodes;
num_elem    = num_nodes -1;
x           = 0:1/(num_nodes-1):1;

input11.coord        = x;
input11.connec(:,1)  = 1:(num_nodes-1);
input11.connec(:,2)  = 2: num_nodes;


%% Global matrix for hydraulic problem


input11.mater(:,1) =     ones(num_elem,1);
input11.mater(:,2) =     ones(num_elem,1);
input11.mater(:,3) = b * ones(num_elem,1);

[Kp,Mp,Cp]  = assembling11(input11);


%% Global matrix for poromechnanical problem
input11.mater(:,1) = (K+4/3*G)*CM *  ones(num_elem,1);
input11.mater(:,2) =        0         *  ones(num_elem,1);
input11.mater(:,3) =        b         *  ones(num_elem,1);

[Ku,~]  = assembling11(input11);

output.Kp = Kp;
output.Ku = Ku;
output.Mp = Mp;
output.Cp = Cp;
end 
