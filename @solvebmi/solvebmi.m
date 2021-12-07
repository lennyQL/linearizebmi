function [gg, vars,outopts] = solvebmi(S, vlist, optg, opts)
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
%     vars:   some data about initial solution and optimal solution
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


%%% setup default option if not determind
% if not opts input
if nargin == 3 || ~isfield(opts,'yalmip')
    %%% set SDP solver
    opt=sdpsettings;
    opt.solver='sedumi';	% 'sedumi' as default SDP solver
    opt.verbose=0;
    opts.yalmip = opt;
end

% select dilation type
if ~isfield(opts,'dilate')
    opts.dilate = false;
end

% debug stdout flag
if ~isfield(opts,'showstep')
    opts.showstep = true;
end

% loop count for option
if ~isfield(opts,'lcmax')
    opts.lcmax = 200;
end


% construct or get Z
isZ = opts.dilate; % checker existence of Z
Z = sdpvar(sizeY(1),sizeY(1));
sizeZ = size(Z);
Z0dummy = sdpvar(sizeZ(1),sizeZ(2));


% debug stdout flag
debug = opts.showstep;

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

% define S size
% class(S)
if isequal(class(S),'cell')
    sizeS = length(S);
elseif isequal(class(S),'string')
    sizeS = 1;
else
    error("varargin{1} must be 'cell' or 'string' class");
end
% sizeS

% S
LMIlist = [];        % constraints
sdpvarnamelist = []; % varnames

for i=1:sizeS
    Fstr = S{i};
    
    if isZ
        [LMIauto,~,gLMI] = linearizebmi(Fstr, {Xstr,Ystr,'Z'}, {'X0dummy','Y0dummy','Z0dummy'});
    else
        [LMIauto,~,gLMI] = linearizebmi(Fstr, vlist, {'X0dummy','Y0dummy'});
    end
    
    
    % LMIauto
    % gLMI.sdpvarname
    
    % declare LMI constraints
    if gLMI.isbmi
        BMImat = LMIauto;
        BMIopt = gLMI;
    else
        LMIlist = [LMIlist, LMIauto<=0];
    end
    
    % sdpvar names in all constraints
    sdpvarnamelist = [sdpvarnamelist, gLMI.sdpvarname];
    
end

% BMImat
% LMIlist
% sdpvarnamelist

sdpvarnamelist = unique(sdpvarnamelist,'stable');



%% search init value by optimize [LMI<t] (t->negative)

%%% epsilon
eps = 1e-6;
upeps = 1e3;
loweps = 1e-3;


%%% desison value 't'
% A dilated LMI constraint is below
% [Q  XN      [tI O
%  GY -G]  -   O  O]  <  O
%

% size of linear term Q in the dilated LMI
sizeQ = size(BMIopt.data.Q);
sizeLMI = size(BMImat);

% decide objective function 't'
t = sdpvar;
tI = t * eye(sizeQ);
alpha = blkdiag(tI, zeros(sizeLMI-sizeQ));

%%% constraint for search init val
LMIinit = [BMImat <= alpha]; % minimize t
% LMIinit = replace(LMIinit,g,1e3);          % make g to big number(constant) 
if isZ
  LMIinit = [LMIinit, Z+Z'>=eps];
%   LMIinit = [LMIinit, Z+Z'>=loweps];
%   LMIinit = [LMIinit, Z+Z'<=upeps];
end


% Append other LMI
LMIinit = [LMIinit, LMIlist];




%%% Init val
% set init val by assign()
% if not assign, set default val
%
% All of these process is just for
% finding max eig of LMIauto as alpha(t)



% Decide others sdpvar default value
firstvaluelist = {};

for i=1:length(sdpvarnamelist)
    name = sdpvarnamelist(i);    % string
    var = evalin('caller',name); % sdpvar
    data = value(var);           % value
    if isnan(data)
        assign(var,zeros(size(data)))
    end
    % Save first init val
    firstvaluelist(i) = {data};
end
Z0first = value(Z);


% Assign intial values
assign(X0dummy,value(X))
assign(Y0dummy,value(Y))
assign(Z0dummy,eye(sizeZ))

% X0, set default value if not exist assign value
X0init = value(X0dummy);
if isnan(X0init)
%     val = zeros(sizeX);
    val = eye(sizeX);
    % val = eye(sizeX) * 1e3;
    assign(X0dummy,val)
%     assign(X,zeros(sizeX))
    X0init = value(X0dummy);
   
end

% Y0
Y0init = value(Y0dummy);
if isnan(Y0init)
    val = zeros(sizeY);
    % val = -ones(sizeY);
    assign(Y0dummy,val)
%     assign(Y,zeros(sizeY))
    Y0init = value(Y0dummy);
end

% Z0
Z0init = value(Z0dummy);
assign(Z,Z0init)

% g
if isnan(value(g))
    assign(g,0);
end




% Assign values to decision values as O
% assign(X,zeros(sizeX))
% assign(Y,zeros(sizeY))
% assign(Z,zeros(sizeZ))


% value(g)
% value(X)
% value(X0dummy)
% value(Y)
% value(Y0dummy)
% value(Z)
% value(Z0dummy)
% value(BMImat)

% Find max eig in LMIauto
eiglmi = eig(value(BMImat));
maxeig = max(eiglmi)
% Decide first value 
% if (>=0)   : search init solution
% elseif (<0): assign val is init solution
tt = maxeig;



%%% stop when t is negative
if debug
  disp("###################################");
  disp('### search for initial solution ###');
  disp("###################################");
  disp("start")
end

ttall=[];   % optimal solutions
lc = 0;
while tt >= 0
  lc=lc+1;
  
  extLMI=LMIinit;
  extLMI=replace(extLMI,X0dummy,X0init);
  extLMI=replace(extLMI,Y0dummy,Y0init);
  extLMI=replace(extLMI,Z0dummy,Z0init);
  

  optimize(extLMI,t,opts.yalmip);

  X0init=double(X);
  Y0init=double(Y);
  Z0init=double(Z);

  
  tt=double(t);
  ttall=[ttall,tt];
  
%   gout=double(g)
  
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
LMI = [BMImat<=-eps*eye(size(BMImat))];
if isZ
  LMI = [LMI, Z+Z'>=eps];
%   LMI = [LMI, Z+Z'>=loweps];
%   LMI = [LMI, Z+Z'<=upeps];
end


% Append other LMI
LMI = [LMI, LMIlist];



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
ggall=[];   % optimal solutions
tmall=[];   % computational times
tStart = tic;

for lc=1:lcmax
  % replace dummy to new optimized val
  extLMI=LMI;
  extLMI=replace(extLMI,X0dummy,X0);
  extLMI=replace(extLMI,Y0dummy,Y0);
  if isZ
    extLMI=replace(extLMI,Z0dummy,Z0);
  end
  
  % optimize by dilated LMI constraits as sufficient conditions of BMI
  optimize(extLMI,g,opts.yalmip);
  tEnd = toc(tStart);
  
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
  
  % show ech step computational time
  tmall=[tmall,tEnd];
  
  if debug
      fprintf('Loop#%03d: %9.4f\n',lc,gg);
  end
  
end


if debug
    disp("end");
end




%% Data of solutions as output
outopts.ggall = ggall; % original objective function
outopts.ttall = ttall; % objective function for initial
outopts.tmall = tmall; % computational time

% % optimal solution
% vars.(inputname(3)) = gg; % original objective value
% vars.(Xstr) = X0;         % decision matrix X
% vars.(Ystr) = Y0;         % decision matrix Y
% if isZ
%   vars.('G') = Z0;       % decision matrix G
% end
% 



% optimal solution
% sdpvarnamelist
for i=1:length(sdpvarnamelist)
    name = sdpvarnamelist(i);
    data = evalin('caller',name);
    vars.(name) = value(data);
end

% vars


% initial solution
vars.([Xstr '0']) = X0init;
vars.([Ystr '0']) = Y0init;
if isZ
  vars.('G0') = Z0init;
end


%% Clear assign value as first intial value

for i=1:length(sdpvarnamelist)
    name = sdpvarnamelist(i);    % string
    var = evalin('caller',name); % sdpvar
    assign(var,firstvaluelist{1,i}); % back to first value
end

assign(Z,Z0first);


