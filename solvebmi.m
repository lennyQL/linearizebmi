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
%                 opt = sdpsetting
%                 opt.solver = 'sedumi'
%                 opts = solvebmiOptions
%                 opts.yalmip = opt
%                 opts.lcmax = 200
%                 opts.method = 1
%                 opts.penalty = 0
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
length(vlist)
try
    Xstr =char(vlist{1})
    Ystr =char(vlist{2})
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
if nargin == 3
    opts = solvebmiOptions;
end


% construct or get Z
switch opts.method % checker existence of Z
    case {0,2}
        isZ = 0;
    case {1,3,4}
        isZ = 1;
    otherwise
        error('not supported method');
end
Z = sdpvar(sizeY(1),sizeY(1),'full');
sizeZ = size(Z);
Z0dummy = sdpvar(sizeZ(1),sizeZ(2),'full');


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

% linearizebmi options
lopts = linearizebmiOptions('method',opts.method);

for i=1:sizeS
    Fstr = S{i};
    
    if isZ
%         Z,is(Z,'scalar')
        [LMIauto,~,gLMI,BMI] = linearizebmi(Fstr, {Xstr,Ystr,'Z'}, {'X0dummy','Y0dummy','Z0dummy'},'',lopts);
    else
        [LMIauto,~,gLMI,BMI] = linearizebmi(Fstr, vlist, {'X0dummy','Y0dummy'},'',lopts);
    end
    
    
    % LMIauto
    % gLMI.sdpvarname
    
    % declare LMI constraints
    if gLMI.isbmi
        BMImat = LMIauto;
        BMIopt = gLMI;
        orgBMI = BMI;
    else
        LMIlist = [LMIlist, LMIauto<=-1e-6];
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
if isZ
  LMIinit = [LMIinit, Z+Z'>=eps];
%   LMIinit = [LMIinit, Z+Z'>=loweps];
%   LMIinit = [LMIinit, Z+Z'<=upeps];
end


% Append other LMI
LMIinit = [LMIinit, LMIlist];

% Limits gamma upperbound
if opts.testg
    LMIinit = [LMIinit, 1e2>=g>=0];
end

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

ttall = tt;   % optimal solutions
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
LMI = [BMImat<=0];
if isZ
  LMI = [LMI, Z+Z'>=eps];
%   LMI = [LMI, Z+Z'>=loweps];
%   LMI = [LMI, Z+Z'<=upeps];
end


% Append other LMI
LMI = [LMI, LMIlist];


% Add regularization term
v = sdpvar(1,1);
pt = opts.penalty;
if ~isequal(pt,0) && ~opts.regterm
    [vxr,vxc]=size(X);
    vx =sdpvar(vxr,vxr,'symmetric');
    if issymmetric(X)
      LMI=[LMI,[vx,triu(X-X0dummy);triu(X-X0dummy)',eye(vxc)]>=0];
    else
      LMI=[LMI,[vx,X-X0dummy;(X-X0dummy)',eye(vxc)]>=0];
    end

    [vyr,vyc]=size(Y);
    vy =sdpvar(vyr,vyr,'symmetric');
    if issymmetric(Y)
      LMI=[LMI,[vy,triu(Y-Y0dummy);triu(Y-Y0dummy)',eye(vyc)]>=0];
    else
      LMI=[LMI,[vy,Y-Y0dummy;(Y-Y0dummy)',eye(vyc)]>=0];
    end

    if isZ
      [vzr,vzc]=size(Z);
      vz =sdpvar(vzr,vzr,'symmetric');
      if issymmetric(Z)
        LMI=[LMI,[vz,triu(Z-Z0dummy);triu(Z-Z0dummy)',eye(vzc)]>=0];
      else
        LMI=[LMI,[vz,Z-Z0dummy;(Z-Z0dummy)',eye(vzc)]>=0];
      end

      LMI=[LMI,v>=trace(vx)+trace(vy)+trace(vz)];
    else
      LMI=[LMI,v>=trace(vx)+trace(vy)];
    end
end


% init val from upper process
X0 = X0init;
Y0 = Y0init;
if isZ
%   Z0 = Z0init;
  Z0 = eye(sizeZ);
end


%%% test
if opts.test
    %%% Initial feasible solutions
    % (K=O is a stabilizing static gain)
    Y0init=zeros(size(Y));
    %%% Calculate initial Lyapunov matrix and H-infinity norm
    initLMI=replace(orgBMI,Y,Y0init);
    optimize([initLMI<=0,LMIlist],g,opts.yalmip);
    X0init=double(X);
    ggsav=double(g);
    %%%
    X0 = X0init;
    Y0 = Y0init;
end


%%
%%% loop
if debug
  disp(" ");
  disp("###################################");
  disp('### search for optimal solution ###');
  disp("###################################");
  disp("start")
end

lcmax=opts.lcmax;
stoptol=opts.stoptol;

vars.('OBJinit') = double(g);
ggsav=double(g)

ggall=ggsav;   % optimal solutions
tmall=0;       % computational times

tStart = tic;

for lc=1:lcmax
  %%% replace dummy to new optimized val
  extLMI=LMI;
  extLMI=replace(extLMI,X0dummy,X0);
  extLMI=replace(extLMI,Y0dummy,Y0);
  if isZ
    extLMI=replace(extLMI,Z0dummy,Z0);
  end
  
  %%% applying regularization term
  if opts.regterm && opts.penalty>0
      lmdc = opts.penalty;
      terms = regterm(lmdc,X,Y,X0,Y0,lc);
      optval = g + terms;
      optimize(extLMI,optval,opts.yalmip);
  else
      optimize(extLMI,g+v*pt,opts.yalmip);
  end
  
  %%% optimize by dilated LMI constraits as sufficient conditions of BMI
%   optimize(extLMI,g+v*pt,opts.yalmip);
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
  
  
  if ggsav-gg<stoptol
      ggsav=gg;
      break;
  end
  ggsav=gg;
  
end


if debug
    disp("end");
end




%% Data of solutions as output
outopts.ggall = ggall; % original objective function
outopts.ttall = ttall; % objective function for initial
% outopts.tmall = [0,tmall]; % computational time
outopts.tmall = tmall;

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




%% Regurilation Terms (Penalty Terms)
% as function
function termval = regterm(lmd,X,Y,X0,Y0,lc)

froX = norm(X-X0,'fro');
froY = norm(Y-Y0,'fro');
maxX0 = max(X0,[],'all');
maxY0 = max(Y0,[],'all');
termX = (lmd * froX^2) / (maxX0 * lc);
termY = (lmd * froY^2) / (maxY0 * lc);
termval = termX + termY;

