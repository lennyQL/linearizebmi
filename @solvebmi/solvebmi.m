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


% get Z if exist
Zstr = '';
if length(vlist) == 3
    Zstr =char(vlist{3});
    Z0str=char(v0list{3});
    Z = evalin('base', Zstr);
    Z0 = evalin('base', Z0str);
    sizeZ = size(Z);
    Z0dummy = sdpvar(sizeZ(1),sizeZ(2));
end
% checker existence of Z
isZ = ~isempty(Zstr);



% linearize bmi
% if nargin == 4
%     % using G
%     LMIauto = linearizebmi(S, vlist, {'X0dummy','Y0dummy'}, G);
% elseif nargin == 3
%     % no G as input
%     LMIauto = linearizebmi(S, vlist, {'X0dummy','Y0dummy'});
% end
if isZ
    LMIauto = linearizebmi(S, vlist, {'X0dummy','Y0dummy','Z0dummy'});
else
    LMIauto = linearizebmi(S, vlist, {'X0dummy','Y0dummy'});
end
% LMI = [LMIauto<=0];
eps = 1e-6;
LMI = [LMIauto<=-eps*eye(size(LMIauto))];
% LMI = [LMIauto<=eps*eye(size(LMIauto))];



%%% run: overbounding approximation method
lcmax=opts.lcmax;	% roop step num
ggall=[];

for lc=1:lcmax
  % replace dummy to new optimized val
  extLMI=LMI;
  extLMI=replace(extLMI,X0dummy,X0);
  extLMI=replace(extLMI,Y0dummy,Y0);
  if isZ
    extLMI=replace(extLMI,Z0dummy,Z0);
  end
  
  % optimize by dilated LMI constraits as sufficient conditions of BMI
  optimize(extLMI,opts.g,opts);

  % update determined val
  X0=double(X);
  Y0=double(Y);
  if isZ
      Z0=double(Z);
  end
  

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
if isZ
  outopt.Z = Z0;
end




