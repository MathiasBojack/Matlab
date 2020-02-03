% This script uses the sequential coupling method for solving the 1D HM
% problem in the following matrix form:
%
%   /                        \  /     \
%   | [M]+ alpha*[Kp] , [C]  |  | Pn+1|
%   |                        |* |     | +
%   | -[C]^T          , [Ku] |  | Un+1|
%   \                        /  \     /