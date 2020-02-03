clear
clc
PoroProperty = PoroElasPara();
P_exact = PoroProperty.gamma*1e6;
% 550256.786500366822
% List_test = log10(2):1:4;
% N_test = round(10.^List_test);% Times discretization
% N_test = [2 10 100 1000];   
% N_test = [1 10 100];
% Iteration = 1:2:5;

N_test = 10;
Iteration = 1;

P =zeros(length(N_test), length(Iteration));
for i = 1: length(N_test)
    
    N = N_test(i);
    for j =1:length(Iteration)
        N_ite = Iteration(j);
        [P(i,j), SV(i,j)] =     MultiDiscretizationFixedStress(N, N_ite);
%         disp([i,j])
    end
end

Error = (abs(P-P_exact)/P_exact);
%%
f1 = figure(1001);
clf
hold on;
f1.Children.YScale ='log';
f1.Children.XScale ='log';
marker = {'+','o','*','.','x','s','d','^','v','>','<','p','h'};
color = {'r','g','b','c','m','k','w'};
for j =1:length(Iteration)
    plt = plot(N_test, Error(:,j));
    txtStart = strcat(strcat( '$10^{',num2str(log10(Error(1  ,j)),1),'}$'));
    txtEnd   = strcat(strcat( '$10^{',num2str(log10(Error(end,j)),1),'}$'));
    text(N_test(1),Error(1,j)*1.1,txtStart,'interpreter','latex')
    text(N_test(end),Error(end,j)*1.1,txtEnd,'interpreter','latex')
    leg1{j} = strcat('Iteration Number =', num2str(Iteration(j)));
    plt.Marker = marker{j};
    plt.Color = color{j};
end
% xlim([1,1e3])
title('Fixed Stress Split','interpreter','latex')
xlabel('Time Discretization Number','interpreter','latex')
ylabel('$\varepsilon = = {\| p-p^{exact}\|}/{\| p^{exact}\|}$',...
    'interpreter','latex')
legend(leg1,'interpreter','latex')
box on; grid on;
% saveas(f1,'ErrorAnalysis_TimeDiscretization_FSS.pdf')

%%
f2 = figure(1002);
clf
hold on;
f2.Children.YScale ='log';

for i =1:length(N_test)
    plt2 = plot(Iteration, Error(i,:));
    leg2{i} = strcat(strcat( '$10^{',num2str(log10(N_test(i)),1),'}$')); 
    plt2.Marker = marker{i+3};
    plt2.Color = color{i+3};
end

legend(leg2,'interpreter','latex')

title('Fixed stress split','interpreter','latex')
xlabel('Iteration Number','interpreter','latex')

ylabel('$\varepsilon = = {\| p-p^{exact}\|}/{\| p^{exact}\|}$',...
    'interpreter','latex')
box on; grid on;
% saveas(f2,'ErrorAnalysis_IterationNumber_FSS.pdf')