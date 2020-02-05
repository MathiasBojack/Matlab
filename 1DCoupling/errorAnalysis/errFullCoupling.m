clear; clc
figure(1)
clf;
xlabel('$X$','interpreter','latex')
ylabel('$P$','interpreter','latex')
hold on
figure(2)
clf
grid on
hold on
xlabel('$T$','interpreter','latex')
ylabel('$U$','interpreter','latex')

commun_tau = 0.1;
%% Input parameters:


sigma0          = -1e6; 
input.L         = 10;
Property        = PoroElasPara();
T0              = Property.mu*Property.CM*input.L^2/Property.k;
input.t         = commun_tau*input.L^2/Property.c;  % physical time

input.T         = input.t/T0; % adimensional time
input.num_nodes = 51;

Property.sigma0 = sigma0;             % analytical solution using dimensional properties
input.sigma0    = sigma0*Property.CM; % adimensional value
input.alpha     = 1;





input.num_tstep = 10;    

input.num_iter = 11;

% full coupling solution, using adimensional parameters
    
[P_fulCoup,U_fulCoup] = adimFullCoupling11(input);
p_fulCoup = P_fulCoup(:,end) / Property.CM;
u_fulCoup      = U_fulCoup(end,:) *input.L;
figure(1)
plot(1-linspace(0,1,input.num_nodes),p_fulCoup/Property.gamma/abs(Property.sigma0),'*')
grid on
figure(2)
plot(linspace(0,input.T, input.num_tstep), u_fulCoup,'*')
    
% sequential coupling solution, using adimensional parameters
   
[P_seqCoup,U_seqCoup] = adimSeqCoupling11(input);

p_seqCoup = P_seqCoup(:,end) / Property.CM;
u_seqCoup      = U_seqCoup(end,:) *input.L;
figure(1)
plot(1-linspace(0,1,input.num_nodes),p_seqCoup/Property.gamma/abs(Property.sigma0),'o')
grid on
figure(2)
plot(linspace(0,input.T, input.num_tstep), u_seqCoup,'o')   



%% analytical solution

% Property.tau   = Property.c * linspace(0,input.t, input.num_tstep)/ input.L^2;
Property.tau = linspace(0,commun_tau, input.num_tstep);
[P_analytc,W_analytc,Z] = Consolidation(Property);
p_analytc = P_analytc(:,end) * Property.gamma*Property.sigma0;
w_analytc = W_analytc * input.L;
figure(1)

hold on;
plot(1-Z,p_analytc/Property.gamma/abs(Property.sigma0),'-')
grid on;

figure(2)
hold on;
plt = plot(linspace(0,input.T, input.num_tstep), w_analytc);   % final value -3.8212e-04 m