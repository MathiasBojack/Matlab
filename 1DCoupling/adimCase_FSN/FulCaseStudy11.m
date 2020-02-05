clear; clc
%% Input parameters:

input.alpha     = 0.6;
input.L = 10;

temp            = PoroElasPara();
T0              = temp.mu*temp.CM*input.L^2/temp.k;
t               = 1e6; 
input.T         = t/T0;  % adimensional time 
input.num_tstep = 10;
input.num_nodes = 101;

X       = 0:1/(input.num_nodes-1):1;
listT   = input.T *(0:1/(input.num_tstep-1):1);

temp.gamma
% vector stress applied at the top X=1, vec_T^D = sigma0 * vec_ex
input.sigma0    = -1e6*temp.CM;


[P,U] = adimFullCoupling11(input);


p = P / temp.CM;
u = U * input.L;


%%

figure(1)
clf
hold on;
for i =2:input.num_tstep
    plot(X, p(:,i)/1e6)
end
ylabel('$P$ [MPa]','interpreter','latex')
xlabel('$X$ [1]','interpreter','latex')
title('Pressure evolution','interpreter','latex')

figure(2)
clf
hold on;
plot(listT, u(input.num_nodes,:))
ylabel('$U$ [m]','interpreter','latex')
xlabel('$T$ [s]','interpreter','latex')
title('Subsidence evolution','interpreter','latex')