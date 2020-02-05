mydir  = pwd;
idcs   = strfind(mydir,'\');
parentDirectory = mydir(1:idcs(end)-1);
 
addpath(parentDirectory)
%%-------------------------------------------------------------------------
%
%                        INPTUT: Poroelastic parameters
%
%   c:              Hydraulic Diffussivity in 1D case       [m2.s^-1]
%   gamma:          Loading coefficient                     [1]
%   L:              Length of the unidimensional sample     [m]
%   t:              Time                                    [s]
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

tList = [ 1e-2  1e0 1e4];

figure(1);
clf
hold on;
    
for i =1:length(tList)
    t = tList(i);
    ConsoPara.sigma0_dot = ConsoPara.sigma0 /t;  
    %%-------------------------------------------------------------------------
    %
    %                        OUTPUT: Pressure and subssidence
    %
    %   P1(z,t):      dimensionless pressure normalized by gamma* sigma0  [1]
    %   U(0,t):       dimensionless Subsidence normalized by L            [1]
    %   tau = c*t/L^2:dimensionneless time                                [1]
    %       Instaneous loading for time t,
    %%-------------------------------------------------------------------------
    tau = PoroProperty.c*t/L^2;
    ConsoPara.tau = tau;
    [P1,W,Z1] = Consolidation(ConsoPara);

    plot(Z1,P1*ConsoPara.gamma*ConsoPara.sigma0/1e6,'k')

    %%-------------------------------------------------------------------------
    %
    %                        OUTPUT: Pressure and subssidence
    %
    %  P2(z,t): dimensionless pressure normalized by L^2*gamma *sigma_dot/2/c[1]
    %   tau = c*t/L:            dimensionneless time                        [1]
    %       Uniform loading rate for time t, dSz/dt = sigma0_dto
    %%-------------------------------------------------------------------------
    %%
    % figure(2)
    % clf
    P_ad = -L^2*ConsoPara.gamma*ConsoPara.sigma0_dot/ConsoPara.c/2;
    [P2, Z2] = UniformLoading(tau);
    plot(1-Z2, P2*P_ad/1e6,'r');
    grid on;
    disp(tau)
end
ylim([0, ceil(ConsoPara.gamma*10)/10])
legText{1} = 'Instaneous Load';
legText{2} = 'Constant Loading Rate';
legend(legText,'interpreter','latex','location','best');

xlabel('Vertical position', 'interpreter','latex')
ylabel('Pressure (MPa)', 'interpreter','latex')
%% 
rmpath(parentDirectory)