function [gg, outopt] = solvebmi(S, vlist, optg, opts)
%SOLVEBMI solve bmi using overbounding approximation method
% this is a wrapper for @linearizebmi
% 
%  input:
%     S:      BMI string
%     vlist:  determinant variable(ex: P,K)
%     optg:   a solution which you want to optimize
%     opts:   some options
%             (ex:
%                 opts = sdpsetting
%                 opts.solver = 'sedumi'
%                 opts.lcmax = 200
%                 opts.dilate = 1
%             )
% 
%  output:
%     gg:     optimal solution
%     outopt: some data about initial solution and optimal solution
% 
% !TODO: 
%   - how to input 'G' and 'opts'
%       - because 'G' will be ignored as input


%% get input value
% input as char
try
    Xstr =char(vlist{1});
    Ystr =char(vlist{2});
%     X0str=char(v0list{1});
%     Y0str=char(v0list{2});
catch 
    error('varargin{2} must be the char list');
end


% determinate value
X = evalin('caller', Xstr);
Y = evalin('caller', Ystr);
% presolve value (init value)
% X0 = evalin('base', X0str);
% Y0 = evalin('base', Y0str);

% get value size
sizeX = size(X);
sizeY = size(Y);

% presolve value(dummy)
X0dummy = sdpvar(sizeX(1),sizeX(2));
Y0dummy = sdpvar(sizeY(1),sizeY(2));


% if not opts input
if nargin == 3
    %%% set SDP solver
    opts=sdpsettings;
    opts.solver='sedumi';	% 'sedumi' as default SDP solver
    opts.verbose=0;
end

% checker existence of Z
if isfield(opts,'dilate') && opts.dilate
    isZ = true;
else
    isZ = false;
end
    
% construct or get Z
Z = sdpvar(sizeY(1),sizeY(1));
sizeZ = size(Z);
Z0dummy = sdpvar(sizeZ(1),sizeZ(2));


% debug stdout flag
debug = ~isfield(opts,'showstep') || opts.showstep;

% loop count for option
if ~isfield(opts,'lcmax')
    opts.lcmax = 200;
end

% target value
g = optg;

%% linearize bmi
% if nargin == 4
%     % using G
%     LMIauto = linearizebmi(S, vlist, {'X0dummy','Y0dummy'}, G);
% elseif nargin == 3
%     % no G as input
%     LMIauto = linearizebmi(S, vlist, {'X0dummy','Y0dummy'});
% end
if isZ
    LMIauto = linearizebmi(S, {Xstr,Ystr,'Z'}, {'X0dummy','Y0dummy','Z0dummy'});
else
    LMIauto = linearizebmi(S, vlist, {'X0dummy','Y0dummy'});
end



%% search init value by optimize [LMI<t] (t->negative)

%%% epsilon
eps = 1e-6;
upeps = 1e3;
loweps = 1e-3;

%%% constraint for search init val
t = sdpvar;
LMIinit = [LMIauto<=t*eye(size(LMIauto))]; % minimize t
LMIinit = replace(LMIinit,g,1e3);          % make g to big number(constant) 
if isZ
  LMIinit = [LMIinit, Z+Z'>=eps];
%   LMIinit = [LMIinit, Z+Z'>=loweps];
%   LMIinit = [LMIinit, Z+Z'<=upeps];
end
LMIinit = [LMIinit, X>=eps];



%%% init val
% X0init = zeros(sizeX);
X0init = eye(sizeX);
% X0init = eye(sizeX) * 1e3;

Y0init = zeros(sizeY);
% Y0init = -ones(sizeY);

Z0init = eye(sizeZ);



%%% stop when t is negative
if debug
  disp("###################################");
  disp('### search for initial solution ###');
  disp("###################################");
  disp("start")
end

tt = 1e3;
lc = 0;
while tt >= 0
  lc=lc+1;
  
  extLMI=LMIinit;
  extLMI=replace(extLMI,X0dummy,X0init);
  extLMI=replace(extLMI,Y0dummy,Y0init);
  extLMI=replace(extLMI,Z0dummy,Z0init);
  

  optimize(extLMI,t,opts);

  X0init=double(X);
  Y0init=double(Y);
  Z0init=double(Z);


  tt=double(t);
  
  % debug output (maybe not necessary)
  if debug
      fprintf('Loop#%03d: %9.4f\n',lc,tt)
  end
  
  
  % loop count upper bound
  if lc == 200
      break
  end
  
end


if debug
    disp("end");
end


% show if initial solution is positive definite
% eig(X0init)



%% run: overbounding approximation method

% constraints
LMI = [LMIauto<=-eps*eye(size(LMIauto))];
if isZ
  LMI = [LMI, Z+Z'>=eps];
%   LMI = [LMI, Z+Z'>=loweps];
%   LMI = [LMI, Z+Z'<=upeps];
end
LMI = [LMI, X>=eps];


% init val from upper proccess
X0 = X0init;
Y0 = Y0init;
if isZ
%   Z0 = Z0init;
  Z0 = eye(sizeZ);
end


%%% loop
if debug
  disp(" ");
  disp("###################################");
  disp('### search for optimal solution ###');
  disp("###################################");
  disp("start")
end

gg = 1e3;
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
  optimize(extLMI,g,opts);

  
  %%% debug
%   if double(g) > gg
%       gg=double(g);
%       ggall=[ggall,gg];
%       break
%   end
  
  % update determined val
  X0=double(X);
  Y0=double(Y);
  if isZ
      Z0=double(Z);
  end
%   eig(X0)
%   eig(Z0)
  
  
  % show each step optimized value
  gg=double(g);
  ggall=[ggall,gg];
  
  if debug
      fprintf('Loop#%03d: %9.4f\n',lc,gg);
  end
  
end


if debug
    disp("end");
end




%% output as options
outopt.ggall = ggall;
outopt.X = X0;
outopt.Y = Y0;
if isZ
  outopt.Z = Z0;
end

outopt.X0 = X0init;
outopt.Y0 = Y0init;
if isZ
  outopt.Z0 = Z0init;
end





