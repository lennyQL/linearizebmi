function opts = linearizebmiOptions(varargin)
%LINEARIZEBMIOPTIONS Create option set for linearizebmi
%
% Syntax
%      opts = linearizebmiOptions
%      opts = linearizebmiOptions('NAME1',VALUE1,'NAME2',VALUE2,...)
%      opts = linearizebmiOptions(OLDOPTIONS,'NAME1',VALUE1,'NAME2',VALUE2,...)


%% Options (default)
%  - append default options here if needed!

%%% Constant factor for decomposition matrix
% Update: G = ^G + dG
%  ^G = t * G_{n-1}
%   H = (1-t) * G_{n-1}
%opts.t = 0.99;
opts.t = 0;

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

