clc
clear
f = fopen('projectPath.txt');
folder = textscan(f,'%s');
fclose(f);
projectPath = fullfile(string(folder{1}));
addpath(projectPath);

%%
fid1 = fopen('RockmatrixResults.dat','rt');

ResultCell = textscan(fid1, '%f %f %f %f %f %f %f', 'Headerlines',1);

Time = ResultCell{1};
Time = Time + Time(3) - Time(2);
Time(1) = 0;

P_Gauss = ResultCell{3};
SX      = ResultCell{4};
SY      = ResultCell{5};
SZ      = ResultCell{6};
Sv      = ResultCell{7};

dSv     = Sv(2:end) - Sv(1:end-1);

dP_Gauss= P_Gauss(2:end) - P_Gauss(1:end-1);


CM_ud = 1.7238e-11; % Ku/M/K
E   = 1e10;
nu = 0.25;
K = E/3/(1-2*nu);
G = E/2/(1+nu);
Kv = K + 4/3*G;
b = 0.833;
B = -b/K/CM_ud; % Skempton

dP_Cal  = B*dSv;

Sv_Cal = K/Kv*SY - 4*G/3/Kv*b*P_Gauss; 




%%
rmpath(projectPath)