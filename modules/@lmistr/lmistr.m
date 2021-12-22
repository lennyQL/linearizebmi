function LMIstr = lmistr(varargin)
%LMISTR 
% linearizebmiの戻り値2LMIstrのオブジェクトクラス
% ex)
%   lmistr(Q,LXN,GYR,G,comp)
%   lmistr(Q,LXN,GYR,G,comp)

% num of inputs
n = length(varargin);
LMIstr.n = n;

% fields of this class
LMIstr.Q = "";
LMIstr.LXN = "";
%
LMIstr.GYR = "";
LMIstr.G = "";
%
LMIstr.O1 = "";
LMIstr.O2 = "";
LMIstr.G1 = "";
LMIstr.G2 = "";
LMIstr.HYR = "";
LMIstr.O3 = "";
LMIstr.H = "";
%
LMIstr.comp = "";

% inputs
if n == 5
    LMIstr.Q = varargin{1};
    LMIstr.LXN = varargin{2};
    LMIstr.GYR = varargin{3};
    LMIstr.G = varargin{4};
    LMIstr.comp = varargin{5};
elseif n == 10
    LMIstr.Q = varargin{1};
    LMIstr.LXN = varargin{2};
    LMIstr.O1 = varargin{3};
    LMIstr.O2 = varargin{4};
    LMIstr.G1 = varargin{5};
    LMIstr.G2 = varargin{6};
    LMIstr.HYR = varargin{7};
    LMIstr.O3 = varargin{8};
    LMIstr.H = varargin{9};
    LMIstr.comp = varargin{10};
else
    error("Number of inputs must be 5 or 10.");
end
    

LMIstr = class(LMIstr, 'lmistr');


end

