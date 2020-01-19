clc
clear
folder = fileread('projectPath.txt');
projectPath = fullfile(folder);
addpath(projectPath);
%%
fid = fopen('injectiondata.txt','rt')

ResultCell = textscan(fid, '%f %f %f %f %f %f');

Time = ResultCell{1};
Pressure = ResultCell{2};
InjectionRate = ResultCell{3};
Disp_West = ResultCell{4};
Disp_North = ResultCell{5};
Disp_Vertical = ResultCell{6};

lwidth = 1.5;

f1 = figure(1);
clf
hold on;
plot(Time,Disp_West-Disp_West(1),'m','linewidth',lwidth)
plot(Time,Disp_North-Disp_North(1),'r','linewidth',lwidth)
plot(Time,Disp_Vertical-Disp_Vertical(1),'g','linewidth',lwidth)


xlabel('Time (s)','interpreter','latex')
ylabel('Displacement ($\mu$m)','interpreter','latex')
box on;
xlim([0, 9000])
pbaspect([1,0.3,1])
legtext{1} = 'West';
legtext{2} = 'North';
legtext{3} = 'Vertical';
l = legend(legtext,'Location','bestoutside');
saveas(f1,'Disp_40.6.tiff')

f2 = figure(2);
clf
plot(Time,InjectionRate)

xlabel('Time (s)','interpreter','latex')
ylabel('InjectionRate (10L/min)','interpreter','latex')
xlim([0, 9000])
pbaspect([1,0.3,1])
box on;

saveas(f2,'InjectionRate_40.6.tiff')


f3 = figure(3);
clf
plot(Time,Pressure)

xlabel('Time (s)','interpreter','latex')
ylabel('Pressure (MPa)','interpreter','latex')
xlim([0, 9000])
pbaspect([1,0.3,1])
box on;

saveas(f3,'Pressure_40.6.tiff')


fclose(fid)
%%
rmpath(projectPath)