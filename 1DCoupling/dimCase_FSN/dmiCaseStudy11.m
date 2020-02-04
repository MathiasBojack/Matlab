clear; clc
%% Input parameters:

input.alpha     = 0.6;
input.L         = 10;
input.T         = 1;
input.num_tstep = 10;
input.num_nodes = 10;

temp    = PoroElasPara();
input.sigma0    = -1e6*temp.CM;

[P,U] = dimFullCoupling(input);