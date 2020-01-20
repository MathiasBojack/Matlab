clear
f = fopen('projectPath.txt');
folder = textscan(f,'%s');
fclose(f);
projectPath = fullfile(string(folder{1}));
addpath(projectPath);


fid1 = fopen('NodalRockmatrixResults.dat','rt');

ResultCell = textscan(fid1, '%f %f %f %f %f %f', 'Headerlines',1);

Time = ResultCell{1};

X = ResultCell{2};
Y = ResultCell{3};
P = ResultCell{4};
Ux = ResultCell{5};
Uy = ResultCell{6};

cutPositionX = 0.5; 
eps = 1e-4;
t = 15;
teps = 1e-5;

plotPoints = (X<(cutPositionX+eps)) & (X>(cutPositionX-eps)) ...
             & (Time<(t+teps)) & (Time>(t-teps));
figure(1)
hold on;
plot(Y(plotPoints), P(plotPoints)/1e6)

figure(2)
hold on;
plot(Y(plotPoints), Uy(plotPoints)/1e6)


rmpath(projectPath)