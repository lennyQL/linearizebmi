function opts = linearizebmiOptions(varargin)
%LINEARIZEBMIOPTIONS Create option set for linearizebmi
%
% Syntax
%      opts = linearizebmiOptions
%      opts = linearizebmiOptions(Name,Value)


%% Options

%%% Constant factor for decomposition matrix
% Update: G = ^G + dG
%  ^G = t * G_{n-1}
%   H = (1-t) * G_{n-1}
opts.t = 0.99;

%%% Method type
% 0: Sebe (2007)
% 1: Sebe (2018)
% 2: Shimomura & Fujii (2005)
% 3: Lee & Hu (2016)
% 4: Ren et al. (2021)
opts.method = 1;


%% Input checks

% function call without argument
if nargin == 0
  return;
end

% error check
if mod(nargin,2)
  error('odd number of arguments');
end

for ac=1:2:nargin
  if ~isfield(opts,varargin{ac})
    error('field "%s" does not exist in options',varargin{ac});
  end
end

% set option values
for i=1:nargin/2
    n = i*2;
    opts.(varargin{n-1}) = varargin{n};
end


end

