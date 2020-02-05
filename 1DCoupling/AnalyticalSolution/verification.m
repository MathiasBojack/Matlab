
%%  Instaneous loading


Property        =  PoroElasPara();
input.t         = 1e4;
input.L         = 10;
Property.tau    = Property.c * input.t/ input.L^2;
Property.sigma0          = -1e6; 


[P,W,Z] = Consolidation(Property)