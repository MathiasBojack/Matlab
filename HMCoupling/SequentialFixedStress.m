% clear
PoroProperty = PoroElasPara();

K   = PoroProperty.K;
b   = PoroProperty.b;
CM  = PoroProperty.CM;
M   = PoroProperty.M;
Ku  = PoroProperty.Ku;
G   = PoroProperty.G;
Kv  = PoroProperty.Kv;
S0  = -1e6; % stress variation [Pa] 

%-------------------------------------------------------------------------%
%                                                                         %
%                           Fixed stress split                            %
%                                                                         %
%-------------------------------------------------------------------------%

N = 22;
Delta_P  = zeros(N,1);
Delta_Ev = zeros(N,1);
Delta_Sv = zeros(N,1);


Delta_P(1)  = 0;
Delta_Ev(1) = 0;
error    = 1e-10;
count    = 0;   
convergence = false;


for i = 1 : N-1
    Delta_Sv(i) = K/Kv*S0 - 4*G/3/Kv*b*Delta_P(i); 
    Delta_P(i+1)  = -M*b/Ku * Delta_Sv(i);
end

Delta_e     = (Delta_P(3:N) - Delta_P(2:N-1))/Delta_P(N); 
Beta_Cal    = Delta_e(2:end)./Delta_e(1:end-1); 
beta        = M*b^2/Ku * 4*G/3/Kv;
Delta_e_Cal = Delta_e(1)*beta.^[1:length(Delta_e)]'; 






f2 = figure(2)
clf
hold on;
plot(1:length(Delta_e), (abs(Delta_e)))
plot(1:length(Delta_e), Delta_e,'*k')
f2.Children.YScale ='log';

f3 = figure(3);
clf
hold on;
p3 = plot(1:N, Delta_P/1e6);
title('Fixed Stress Split','interpreter','latex')
xlabel('Iteration Number','interpreter','latex')
ylabel('Terzaghi Pressure (MPa)','interpreter','latex')
text(0,Delta_P(end)/1e6, num2str(Delta_P(end)/1e6))
text(1,Delta_P(2)/1e6, num2str(Delta_P(2)/1e6))



plot([0,N],[1,1]*Delta_P(N)/1e6,'--')
grid on;
saveas(f3,'TerzaghiFixedStress.tiff')

