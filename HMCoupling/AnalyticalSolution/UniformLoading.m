%%-------------------------------------------------------------------------
%
%         This document calculates the pressure evolution of the Terzaghi
%         Problem. Solution comes from wang,2000. Page 137.
%
%%-------------------------------------------------------------------------

%%-------------------------------------------------------------------------
%
%                        INPTUT: Poroelastic parameters
%
%   c:              Hydraulic Diffussivity in 1D case       [m2.s^-1]
%   kv = K + 4/3*G Vertical incompressibility               [Pa]
%   b:              Coefficient de Biot                     [1]
%   gamma:          Loading coefficient                     [1]
%   tau = c*t/L^2:            dimensionneless time          [1]
%   L:              Length of the unidimensional sample     [m]
%   t:              Time 
%   sigma0_dot:     Applied stress at the top of the sample [Pa.s^-1]        
%
%           dp/dt 
%%-------------------------------------------------------------------------



%%-------------------------------------------------------------------------
%
%                        OUTPUT: Pressure and subssidence
%
%   Dimensionless pressure: L^2*gamma*sigma0_dot/c/2
%   P(z,tau)= p(z,tau)/gamma/sigma0_dot/L^2:	dimensionless pressure  [1]
%   U(0,tau)= u(0,tau)/L:               Subsidence              [1]
%
%%-------------------------------------------------------------------------

function [P, Z] = UniformLoading(tau)
% tau     = ConsoPara.tau;
% b       = ConsoPara.b;
% Kv      = ConsoPara.Kv;
% gamma   = ConsoPara.gamma;
% sigma0_dot  = ConsoPara.sigma0_dot; %  sigma0<0, means increasing compression.


Z = [0:0.001:1]';
P = (1 - (1-Z).^2);
% P= zeros(length(Z),1);
for m = 0:1:100
    P = P - 32/pi^3* (-1)^m/(2*m+1)^3 * exp(-(2*m+1)^2*pi^2/4*tau) * cos((2*m+1)*pi/2*(1-Z));
    % cm = b/Kv;
end

