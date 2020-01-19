clc
clear
folder = fileread('projectPath.txt');
projectPath = fullfile(folder);
addpath(projectPath);
%%

%-------------------------------------------------------------------------%
%                                                                         %
%                          1  Experiment BFS2 40.6m                       %
%                                                                         %
%-------------------------------------------------------------------------%

fid1 = fopen('jgrb52263-sup-0003-2017jb014384-ts02.txt','rt');

lwidth = 1.5;

% time(s),Pressure(MPa),displacement(north), displacement(west), ...
% displacement(vertical)
ResultCell1 = textscan(fid1, '%f %f %f %f %f', 'Headerlines',1,...
    'Delimiter',',');

Time1        = ResultCell1{1};        % [s]
Pressure1    = ResultCell1{2};        % [bar] = 0.1[MPa]
Disp_North1  = ResultCell1{3};        % [1e-6m]
Disp_West1   = ResultCell1{4};
Disp_Ver1    = ResultCell1{5};
Disp_Norm1   = vecnorm([Disp_North1, Disp_West1, Disp_Ver1],2,2);
%-------------------------------------------------------------------------%
%                                                                         %
%                          1.1  Pressure                                  %
%                                                                         %
%-------------------------------------------------------------------------%

f11 = figure(11);
clf;
plot(Time1, Pressure1/10,'linewidth',lwidth)

xlabel('Time (s)','interpreter','latex')
ylabel('Pressure (MPa)','interpreter','latex')
xlim([0, 9000])
pbaspect([1,0.3,1])
box on;

saveas(f11,'Pressure_40.6.tiff')
%-------------------------------------------------------------------------%
%                                                                         %
%                          1.2  Displacement                              %
%                                                                         %
%-------------------------------------------------------------------------%
f12 = figure(12);
clf;
hold on;
plot(Time1, Disp_North1,'r','linewidth',lwidth)
plot(Time1, Disp_West1,'m','linewidth',lwidth)
plot(Time1, Disp_Ver1,'g','linewidth',lwidth)


xlabel('Time (s)','interpreter','latex')
ylabel('Displacement ($\mu$m)','interpreter','latex')
box on;
xlim([0, 9000])
pbaspect([1,0.3,1])
legtext{1} = 'North';
legtext{2} = 'West';
legtext{3} = 'Vertical';
l12 = legend(legtext,'Location','bestoutside');
saveas(f12,'Disp_40.6.tiff')
%-------------------------------------------------------------------------%
%                                                                         %
%                          1.3 Normal Displacement                        %
%                                                                         %
%-------------------------------------------------------------------------%
f13 = figure(13);
clf;
plot(Time1, Disp_Norm1,'r','linewidth',lwidth)
xlabel('Time (s)','interpreter','latex')
ylabel('Displacement ($\mu$m)','interpreter','latex')
box on;
xlim([0, 9000])
pbaspect([1,0.3,1])
l13 = legend('Norm Displacement');
saveas(f13,'NormDisp_40.6.tiff')
fclose(fid1);
%% Experiment BFS2 37.2m
%-------------------------------------------------------------------------%
%                                                                         %
%                          2  Experiment BFS2 37.2m                       %
%                                                                         %
%-------------------------------------------------------------------------%

fid2 = fopen('jgrb52263-sup-0002-2017jb014384-ts01.txt','rt');
ResultCell2 = textscan(fid2, '%f %f %f %f', 'Headerlines',1,...
    'Delimiter',',');

% time(s),Pressure(bar),deofrmation (north), deformation(west)
Time2        = ResultCell2{1};        % [s]
Pressure2    = ResultCell2{2};        % [bar] = 0.1[MPa]
Def_North2  = ResultCell2{3};        % 
Def_West2   = ResultCell2{4};





f21 = figure(21)
clf
plot(Time2, Pressure2/10,'linewidth',lwidth)
xlabel('Time (s)','interpreter','latex')
ylabel('Pressure (MPa)','interpreter','latex')
title('Pressure evolution - 37.2m','interpreter','latex')
box on;
pbaspect([1,0.3,1])
ylim([0, 3])
saveas(f21,'Pressure_37.2.tiff')


f22 = figure(22);
clf;
hold on;
plot(Time2,Def_North2,'r','linewidth',lwidth)
plot(Time2,Def_West2, 'm','linewidth',lwidth)
ylim([-3, 2]*1e-3)


xlabel('Time (s)','interpreter','latex')
ylabel('$\varepsilon$','interpreter','latex')
title('Deformation evolution - 37.2m','interpreter','latex')

box on;
% xlim([0, 800])
pbaspect([1,0.3,1])
legtext2{1} = 'North';
legtext2{2} = 'West';
l12 = legend(legtext2,'Location','bestoutside');
saveas(f22,'Disp_37.2.tiff')

%% Experiment BFS2 44.65m
fid3 = fopen('jgrb52263-sup-0004-2017jb014384-ts03.txt','rt');
ResultCell3 = textscan(fid3, '%f %f %f %f', 'Headerlines',1,...
    'Delimiter',',');

% time(s),Pressure(bar),deofrmation (north), deformation(west)
Time3        = ResultCell3{1};        % [s]
Pressure3    = ResultCell3{2};        % [bar] = 0.1[MPa]
Def_North3  = ResultCell3{3};        % 
Def_West3   = ResultCell3{4};





f31 = figure(31)
clf
plot(Time3, Pressure3/10,'linewidth',lwidth)
xlabel('Time (s)','interpreter','latex')
ylabel('Pressure (MPa)','interpreter','latex')
title('Pressure evolution - 44.65m','interpreter','latex')
box on;
% xlim([0, 800])
pbaspect([1,0.3,1])
ylim([0, 3])

saveas(f31,'Pressure_44.65.tiff')


f32 = figure(32);
clf;
hold on;
plot(Time3,Def_North3,'r','linewidth',lwidth)
plot(Time3,Def_West3, 'm','linewidth',lwidth)
ylim([-3, 2]*1e-3)


xlabel('Time (s)','interpreter','latex')
ylabel('$\varepsilon$','interpreter','latex')
title('Deformation evolution - 44.65m','interpreter','latex')

box on;
% xlim([0, 800])
pbaspect([1,0.3,1])
legtext2{1} = 'North';
legtext2{2} = 'West';
l13 = legend(legtext2,'Location','bestoutside');

saveas(f32,'Disp_44.65.tiff')

%% Experiment BFS1 47.2m
fid4 = fopen('jgrb52263-sup-0005-2017jb014384-ts04.txt','rt');
ResultCell4 = textscan(fid4, '%f %f %f %f', 'Headerlines',1,...
    'Delimiter',',');

% time(s),Pressure(bar),deofrmation (north), deformation(west)
Time4        = ResultCell4{1};        % [s]
Pressure4    = ResultCell4{2};        % [bar] = 0.1[MPa]
Def_North4  = ResultCell4{3};        % 
Def_West4   = ResultCell4{4};





f41 = figure(41)
clf
plot(Time4, Pressure4/10,'linewidth',lwidth)
xlabel('Time (s)','interpreter','latex')
ylabel('Pressure (MPa)','interpreter','latex')
title('Pressure evolution - 47.2m','interpreter','latex')
box on;
% xlim([0, 800])
pbaspect([1,0.3,1])
ylim([0, 3])
saveas(f41,'Pressure_47.2.tiff')


f42 = figure(42);
clf;
hold on;
plot(Time4,Def_North4,'r','linewidth',lwidth)
plot(Time4,Def_West4, 'm','linewidth',lwidth)


xlabel('Time (s)','interpreter','latex')
ylabel('$\varepsilon$','interpreter','latex')
title('Deformation evolution - 47.2m','interpreter','latex')

box on;
% xlim([0, 800])
pbaspect([1,0.3,1])
legtext2{1} = 'North';
legtext2{2} = 'West';
l13 = legend(legtext2,'Location','bestoutside');
% ylim([-3, 2]*1e-3)

saveas(f42,'Disp_47.2.tiff')

%%
rmpath(projectPath)