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
ConsoPara.sigma0 = -1e6;
L = 10;
t = 1e-1;
ConsoPara.sigma0_dot = ConsoPara.sigma0 /t;  
%%-------------------------------------------------------------------------
%
%                        OUTPUT: Pressure and subssidence
%
%   P1(z,t):      dimensionless pressure normalized by gamma* sigma0  [1]
%   U(0,t):       dimensionless Subsidence normalized by L            [1]
%   tau = c*t/L:            dimensionneless time                      [1]
%       Instaneous loading for time t,
%%-------------------------------------------------------------------------
tau = PoroProperty.c*t/L;
ConsoPara.tau = tau;
[P1,W,Z1] = Consolidation(ConsoPara);

figure(1);
clf
hold on;
hold on;
plot(Z1,P1*ConsoPara.gamma*ConsoPara.sigma0/1e6,'k*')

%%-------------------------------------------------------------------------
%
%                        OUTPUT: Pressure and subssidence
%
%   P2(z,t):     dimensionless pressure normalized by L^2*gamma *sigma_dot/c[1]
%   tau = c*t/L:            dimensionneless time                        [1]
%       Uniform loading rate for time t, dSz/dt = sigma0_dto
%%-------------------------------------------------------------------------
figure(2)
clf

[P2, Z2] = UniformLoading(ConsoPara);
plot(Z2, P2,'ro')L^2*ConsoPara.sigma0_dot * ConsoPara.gamma/ConsoPara.c /1e6