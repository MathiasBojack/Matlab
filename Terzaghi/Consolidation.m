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
%   kv = K + 4/3*G Vertical incompressibility              [Pa]
%   b:              Coefficient de Biot                     [1]
%   gamma:          Loading coefficient                     [1]
%   tau = c*t/L^2:            dimensionneless time            [1]
%   L:              Length of the unidimensional sample     [m]
%   t:              Time 
%   sigma0:         Applied stress at the top of the sample         
%
%%-------------------------------------------------------------------------



%%-------------------------------------------------------------------------
%
%                        OUTPUT: Pressure and subssidence
%
%   P(z,tau)= p(z,tau)/gamma/sigma0:	dimensionless pressure  [1]
%   U(0,tau)= u(0,tau)/L:               Subsidence              [1]
%
%%-------------------------------------------------------------------------

function [P,W,Z] = Consolidation(ConsoPara)
tau     = ConsoPara.tau;
b       = ConsoPara.b;
Kv      = ConsoPara.Kv;
gamma   = ConsoPara.gamma;
sigma0  = ConsoPara.sigma0;


Z = [0:0.01:1]';
P = zeros(length(Z),1);
W = b/Kv * gamma*sigma0;
for m = 0:1:100
    P = P + 4/pi / (2*m+1) * exp(-(2*m+1)^2*pi^2*tau/4) * sin((2*m+1)*pi/2*Z);
    W = W - b/Kv*gamma*sigma0 * 8/pi^2 / (2*m+1)^2 * exp(-(2*m+1)^2*pi^2*tau/4);
end

