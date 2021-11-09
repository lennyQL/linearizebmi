function [gg, outopt] = solvebmi(S, vlist, v0list, opts)
%SOLVEBMI solve bmi using overbounding approximation method
% this is a wrapper for @linearizebmi
% !TODO: 
%   - how to input 'G' and 'opts'
%       - because 'G' can ignore as input


% get input value
% input as char
Xstr =char(vlist{1});
Ystr =char(vlist{2});
X0str=char(v0list{1});
Y0str=char(v0list{2});

% determinate value
X = evalin('base', Xstr);
Y = evalin('base', Ystr);
% presolve value (init value)
X0 = evalin('base', X0str);
Y0 = evalin('base', Y0str);

% get value size
sizeX = size(X);
sizeY = size(Y);

% presolve value(dummy)
X0dummy = sdpvar(sizeX(1),sizeX(2));
Y0dummy = sdpvar(sizeY(1),sizeY(2));



% linearize bmi
% if nargin == 4
%     % using G
%     LMIauto = linearizebmi(S, vlist, {'X0dummy','Y0dummy'}, G);
% elseif nargin == 3
%     % no G as input
%     LMIauto = linearizebmi(S, vlist, {'X0dummy','Y0dummy'});
% end
LMIauto = linearizebmi(S, vlist, {'X0dummy','Y0dummy'});
LMI = [LMIauto<=0];



%%% run: overbounding approximation method
lcmax=opts.lcmax;	% roop step num
ggall=[];

for lc=1:lcmax
  extLMI=LMI;
  extLMI=replace(extLMI,X0dummy,X0);
  extLMI=replace(extLMI,Y0dummy,Y0);

  optimize(extLMI,opts.g,opts);

  X0=double(X);
  Y0=double(Y);

  % show each step optimized value
  gg=double(opts.g);
  ggall=[ggall,gg];
  
  if ~isfield(opts,'showstep') || opts.showstep
      fprintf('Loop#%03d: %9.4f\n',lc,gg);
  end
  
end


outopt.ggall = ggall;
outopt.X = X0;
outopt.Y = Y0;




