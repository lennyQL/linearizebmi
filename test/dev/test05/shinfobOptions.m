% shinfobOptions creates option set for hinfob1 and hinfob2
%
% Syntax
%      opts = shinfobOptions()
%      opts = shinfobOptions(Name,Value)

function opts = shinfobOptions(varargin)

%%% default option set
%% flags
% flags for control
opts.flag_dual=0;
opts.flag_fix_h=0;
opts.flag_speed=0;
opts.flag_analysis=0;		% norm analysis in main loop
opts.lcmax=1e3;			% maximum number of iteration
opts.max_retry=20;		% maximum number of retry
opts.flag_kpnt=0;		% penalty for structural constraint
opts.kpnt_msk=[];		% mask for penalty on controller parameters
opts.kpnt_pnt=0;		% penalty on controller parameters
opts.kmsk=[];			% structural constraint on controller
opts.flag_k_stable=0;		% design a stable controller
% flags for numerical control
opts.flag_pbound_init=1;
opts.flag_plbound_init=1;
opts.flag_pbound=1;
opts.flag_plbound=1;


%% numerical options
% main options
opts.gmin=0;			% target H-infinity norm
opts.ginit=1;			% initial guess of H-infinity norm
opts.stop=1e-6;			% stopping threshold
opts.init_sc.in=[];		% system input-side initial scaling
opts.init_sc.out=[];            % system output-side initial scaling

% tolerances and bounds
%opts.macheps=eps;
%opts.tolsing=sqrt(macheps);
%opts.toleig=macheps^(2/3);

opts.feasrad=1e10;
%opts.lmi_options=[1e-6,1e3,opts.feasrad,20,1];
opts.lmi_options=[1e-6,1e3,opts.feasrad,40,1];
opts.lmieps=1e-6;
opts.pbound=1e10;

opts.pnt_p_init=1e-4;	% add penalty trace(P)*pnt_p to objective function
opts.pnt_p_init=0;
opts.pnt_p=0;

opts.k_stable_pnt=1e-6;		% penalty for controller Lyapunov matrix
opts.k_stable_plbound=1;	% lower bound for controller Lypunov matrix

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
    error('filed %s does not exist in options',varargin{ac});
  end
end

% set option values
opts=setfield(opts,vargin);

end

