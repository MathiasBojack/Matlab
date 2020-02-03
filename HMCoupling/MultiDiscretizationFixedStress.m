function [Pressure, SV] = MultiDiscretizationFixedStress(N, N_Iteration)
%%
PoroProperty = PoroElasPara();

K   = PoroProperty.K;
b   = PoroProperty.b;
CM  = PoroProperty.CM;
M   = PoroProperty.M;
Ku  = PoroProperty.Ku;
G   = PoroProperty.G;
Kv  = PoroProperty.Kv;
S0  = -1e6; % stress variation [Pa] 


%% 

% There is only one iteration from the hydraulic to the mechanical Problem
% Using the fixed stress method
%
% N = 10000; %% Loading step
% N_Iteration = 1;





% With the definition of intial values, dSv, dSz, dP are extended to N_Step
% length.

Sz_initial = 0;
Sv_initial = 0;
P_initilal  = 0;

dSz = ones(N,1)*(S0-Sz_initial)/N;   % dSz(k-1) = Sz(k) - Sz(k-1)
dSv = zeros(N,1);                    % dSv(k-1) = Sv(k) - Sv(k-1)
dP  = zeros(N,1);                    % dP(k-1) = P(k) - P(k-1)
ScrH = zeros(N,1);
dSvk_1_i = zeros(N_Iteration,1);
dPk_1_i  = zeros(N_Iteration,1);

P = zeros(N+1,1);
P(1) = P_initilal;
Sv = zeros(N+1,1);
Sv(1) = Sv_initial;

for k = 2 : N+1
    
    dPk_1_i(1) = dP(k-1);   
     
    for i = 1:N_Iteration-1
        % Fixed stress method
        dSvk_1_i(i)  = K/Kv*dSz(k-1) - 4*G/3/Kv*b*dPk_1_i(i);
        dPk_1_i(i+1) = - M*b/Ku * dSvk_1_i(i);       
    end
    dSvk_1_i(N_Iteration)  = K/Kv*dSz(k-1) - 4*G/3/Kv*b*dPk_1_i(N_Iteration);
    % Pass the final result of in-step iteration to the next time step  
    dSv(k-1) = dSvk_1_i(end); 
%     dP(k-1) = - M*b/Ku * dSv(k-1);
% for comparison with DisRoc
    dP(k) = - M*b/Ku * dSv(k-1);
    P(k) = P(k-1) + dP(k-1);
    Sv(k) = Sv(k-1) + dSv(k-1);
    ScrH(k-1) = -dSv(k-1)*b/K*N; % dt = 1/N
end

Pressure = P(end);
SV = Sv(end);
end
% 
% figure(101)
% plot(1:N+1, P/abs(S0))
% text(N+1,P(end)/abs(S0),num2str(P(end)/abs(S0)))