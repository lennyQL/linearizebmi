function options = solvebmisettings(varargin)
%SOLVEBMISETTINGS Create solvebmi option structure (default)

options.yalmip = sdpsettings;
% options.sdpsettings = sdpsettings;

options.lcmax = 2e2;

options.stoptol = 5e-7;

options.showstep = 1;

% dilation types
% - false: decomposition matrix is constant (2x2 blocks LMI)
% - true:  decomposition matrix is decision value (3x3 blocks LMI)
options.dilate = 0;

% regularetion term
% options.regterm = 0;
% Penalty term factor:
% - this is not bool; but numeric
% - if 0, there is no penaty terms
options.penalty = 0;

end

