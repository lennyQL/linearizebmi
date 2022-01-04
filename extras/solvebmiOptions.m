function opts = solvebmiOptions(varargin)
%SOLVEBMIOPTIONS Create option set for solvebmi
%
% Syntax
%      opts = solvebmiOptions
%      opts = solvebmiOptions('NAME1',VALUE1,'NAME2',VALUE2,...)
%      opts = solvebmiOptions(OLDOPTIONS,'NAME1',VALUE1,'NAME2',VALUE2,...)

%% Options (default)
%  - append default options here if needed!

opts.yalmip = sdpsettings;
% options.sdpsettings = sdpsettings;

opts.lcmax = 2e2;

opts.stoptol = 5e-7;

opts.showstep = 1;

%%% dilation types
% - false: decomposition matrix is constant (2x2 blocks LMI)
% - true:  decomposition matrix is decision value (3x3 blocks LMI)
opts.dilate = 0;

%%% regularetion term
opts.regterm = 0;
%%% Penalty term factor:
% - this is not bool; but numeric
% - if 0, there is no penaty terms
opts.penalty = 0;



%% Input checks

% function call without argument
if nargin == 0
  return;
end

% input pre options
namestart = 0;
if isa(varargin{1},'struct')
    preopts = varargin{1};
    namestart = 1;
    % how to append pre input options?
    namelist = fieldnames(preopts);
    for i=1:length(namelist)
        opts.(namelist{i}) = preopts.(namelist{i});
    end    
end

% error check
if mod(nargin-namestart,2)
  error('odd number of arguments');
end

for ac=1+namestart:2:nargin
  if ~isfield(opts,varargin{ac})
    error('field "%s" does not exist in options',varargin{ac});
  end
end


% set option values
for i=1:(nargin-namestart)/2
    n = i*2+namestart;
    opts.(varargin{n-1}) = varargin{n};
end


end

