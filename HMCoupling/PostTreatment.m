%%-------------------------------------------------------------------------
%
%                        INPTUT: Poroelastic parameters
%
%   c:              Hydraulic Diffussivity in 1D case       [m2.s^-1]
%   gamma:          Loading coefficient                     [1]
%   L:              Length of the unidimensional sample     [m]
%   t:              Time 
%   sigma0:         Applied stress at the top of the sample         
%
%%-------------------------------------------------------------------------

PoroProperty    = PoroElasPara();
ConsoPara.c     = PoroProperty.c;
ConsoPara.gamma = PoroProperty.gamma;
ConsoPara.b     = PoroProperty.b;
ConsoPara.Kv    = PoroProperty.Kv;
ConsoPara.sigma0 = 1e6;
L = 10;
t = 10;
%%-------------------------------------------------------------------------
%
%                        OUTPUT: Pressure and subssidence
%
%   P(z,t):       dimensionless pressure                    [Pa]
%   U(0,t):       dimensionless Subsidence                  [m]
%   tau = c*t/L:            dimensionneless time            [1]
%
%%-------------------------------------------------------------------------
tau = PoroProperty.c*t/L;
ConsoPara.tau = tau;
[P,W,Z] = Consolidation(ConsoPara);

figure();
hold on;
plot(Z,P*ConsoPara.gamma)