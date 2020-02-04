
clear; clc
%% Input parameters:


sigma0          = -1e6; 
Property        =  PoroElasPara();

input.T         = 1e1;
input.num_nodes = 101;
input.sigma0    = sigma0*Property.CM; % adimensional value
input.alpha     = 0.65;
input.L         = 10;


Property.sigma0 = sigma0;

list_tstep = [ 1e1 1e2 1e3];
list_niter = [ 1:2:10];

leapfrog_err_fulCoup = zeros(length(list_tstep),1);


%% analytical solution

Property.tau   = Property.c * input.T / input.L^2;
[P_analytc,W_analytc,Z] = Consolidation(Property);
p_analytc = P_analytc * Property.gamma*Property.sigma0;

figure(1)
clf;
hold on;
plot(Z,p_analytc,'r')

%% full coupling solution
for i = 1: length(list_tstep)
    
    input.num_tstep = list_tstep(i);    
    [P_fulCoup,U_fulCoup] = adimFullCoupling11(input);
    p_fulCoup(:,i) = P_fulCoup(:,end) / Property.CM;
    
    figure(1)
    plot(linspace(0,1,input.num_nodes),p_fulCoup(:,i),'k')
    

%     plot
    
%     leapfrog_err_fulCoup = ;
        
%     for k = 1: length(list_niter) 
%         
%         input.num_iter = list_niter(j) ;
%         
% 
%         
%     end
%     leapfrog_err_fulCoup(:,i) = norm(())/norm()
end
