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

dSz = ones(N,1)*(S0-Sz_initial)/N;
dSv = zeros(N,1);
dP  = zeros(N,1);

dSv_k_1 = zeros(N_Iteration,1);
dP_k_1  = zeros(N_Iteration,1);

P = zeros(N+1,1);
P(1) = P_initilal;
Sv = zeros(N+1,1);
Sv(1) = Sv_initial;

for k = 2 : N+1
    
    dP_k_1(1) = dP(k-1);   
%     
%     if N_Iteration == 1
%         dSv_k(1) = K/Kv*dSz(k) - 4*G/3/Kv*b*dP_k(1);
%     end
     
    for i = 1:N_Iteration-1
        % Fixed stress method
        dSv_k_1(i)  = K/Kv*dSz(k-1) - 4*G/3/Kv*b*dP_k_1(i);
        dP_k_1(i+1) = - M*b/Ku * dSv_k_1(i);       
    end
    dSv_k_1(N_Iteration)  = K/Kv*dSz(k-1) - 4*G/3/Kv*b*dP_k_1(N_Iteration);
    % Pass the final result of in-step iteration to the next time step  
    dSv(k-1) = dSv_k_1(end); 
    dP(k-1) = - M*b/Ku * dSv(k-1);
    P(k) = P(k-1) + dP(k-1);
    Sv(k) = Sv(k-1) + dSv(k-1);
end

Pressure = P(end);
SV = Sv(end);
end
% 
% figure(101)
% plot(1:N+1, P/abs(S0))
% text(N+1,P(end)/abs(S0),num2str(P(end)/abs(S0)))