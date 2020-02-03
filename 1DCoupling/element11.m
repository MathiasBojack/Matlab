function [K,M,C] = element11(L)

% 1D bar element for 2nd order PDE in x1
% Linear element using the following shape function:
% 
%       N1(x1) = 1-x1/L; N2(x1) = x1/L
% 
% where : L is the length of the element.

% The interpolation for the unknown variable in the interval [0,L] could be
% written in the following matrix form:
%
%       Uh(x_1) = [N]*{U^e}
%
% while its graident could be written as:
%
%       grad_Uh(x_1) = [B]*{U^e}
%
% where the gradient matrix is defined by:
%
%       [B] =1/L*[-1, 1]
%
% The elemententary matrix K,M,C are defined by:
%
%       K = int_0^L{[B]^T*[B]*dx1}
%       M = int_0^L{[N]^T*[N]*dx1}
%       C = int_0^L{[N]^T*[B]*dx1}

% Element subroutine to construct local stiffness matrix from local data

 
M = L/6*[ 2  1; 1,2];
K = 1/L*[ 1 -1;-1,1];
C = 1/2*[-1  1;-1,1];

end