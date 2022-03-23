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
%   - how to input constant 'G' data
%       - because 'G' will be ignored as input
%       - in some method, fixed 'G' is required


%% get input value
% input as char

% define S size (length of input constraints)
% class(S)
if isequal(class(S),'cell')
    sizeS = length(S);
elseif isequal(class(S),'string')
    sizeS = 1;
else
    error("varargin{1} must be 'cell' or 'string' class");
end
% sizeS


% (vlen): Number of decision variables set (Number of input BMIs)
vlen = 0;
for i=1:length(vlist)
    if length(vlist{i}) == 1
        vlen = 1;
        sizeV = 1;
        break
    else
        sizeV = length(vlist);
        if ~isempty(vlist{i})
            vlen = vlen + 1;
        end
    end
end


% Size S and vlist must be same
if ~isequal(sizeS,sizeV) && sizeV > 1
    error("Length of varargin{1} and varargin{2} must be the same")
end


if vlen == 1 && isa(vlist{1},'char')
    vlist = {vlist};
end

% vlist for bmi
bmivlist = {};
for i=1:sizeV
    if ~isempty(vlist{i})
        bmivlist = cat(1,bmivlist,vlist(i));
    end
end


% (vstruct): Info about decision variables set
% !TODO: vstruct need to be a class (maybe)
vstruct(1:vlen) = struct;
for i=1:vlen
    % var str
    try
        Xstr = char(bmivlist{i}{1});
        Ystr = char(bmivlist{i}{2});
    catch 
        error('varargin{2} must be the char list (length:2)');
    end
    vstruct(i).Xstr = Xstr;
    vstruct(i).Ystr = Ystr;
    % determinate value data
    X = evalin('caller', Xstr);
    Y = evalin('caller', Ystr);
    vstruct(i).X = X;
    vstruct(i).Y = Y;
    % get value size
    sizeX = size(X);
    sizeY = size(Y);
    vstruct(i).sizeX = sizeX;
    vstruct(i).sizeY = sizeY;
    % presolve value (dummy)
    if issymmetric(X)
        X0dummy = sdpvar(sizeX(1),sizeX(2),'symmetric');
    else
        X0dummy = sdpvar(sizeX(1),sizeX(2),'full');
    end
    if issymmetric(Y)
        Y0dummy = sdpvar(sizeY(1),sizeY(2),'symmetric');
    else
        Y0dummy = sdpvar(sizeY(1),sizeY(2),'full');
    end
    vstruct(i).X0dummy = X0dummy;
    vstruct(i).Y0dummy = Y0dummy;
    % info about Z
    if opts.method == 3
        Z = sdpvar(sizeY(1),sizeY(1),'symmetric');
        sizeZ = size(Z);
        Z0dummy = sdpvar(sizeZ(1),sizeZ(2),'symmetric');
    else
        Z = sdpvar(sizeY(1),sizeY(1),'full');
        sizeZ = size(Z);
        Z0dummy = sdpvar(sizeZ(1),sizeZ(2),'full');
    end
    vstruct(i).Z = Z;
    vstruct(i).sizeZ = sizeZ;
    vstruct(i).Z0dummy = Z0dummy;
    % info about M
    M = sdpvar(sizeY(1),sizeY(1),'symmetric');
    sizeM = size(M);
    M0dummy = sdpvar(sizeM(1),sizeM(2),'symmetric');
    vstruct(i).M = M;
    vstruct(i).sizeM = sizeM;
    vstruct(i).M0dummy = M0dummy;
    
end
% vstruct(1)
% vstruct(2)


%%% setup default option if not determind
% if not opts input
if nargin == 3
    opts = solvebmiOptions;
end


% construct or get Z
switch opts.method % checker existence of Z
    case {0,2}
        isZ = 0;
        isM = 0;
    case {1,3}
        isZ = 1;
        isM = 0;
    case 4
        isZ = 1;
        isM = 1;
    otherwise
        error('not supported method');
end


% debug stdout flag
debug = opts.showstep;

% target value
g = optg;

%% linearize bmi

% S
LMIlist = [];        % constraints
sdpvarnamelist = []; % varnames

% linearizebmi options
lopts = linearizebmiOptions('method',opts.method,'t',opts.t);

vstep = 1;  % step of vstruct(corresponding to BMI)
bminum = 0; % number of BMIs
for i=1:sizeS
    Fstr = S{i};
    
    Xstr = vstruct(vstep).Xstr;
    Ystr = vstruct(vstep).Ystr;
    Z = vstruct(vstep).Z;
    M = vstruct(vstep).M;
    X0dummy = vstruct(vstep).X0dummy;
    Y0dummy = vstruct(vstep).Y0dummy;
    Z0dummy = vstruct(vstep).Z0dummy;
    M0dummy = vstruct(vstep).M0dummy;
    
    
    if isM
        try
            [LMIauto,~,gLMI,BMI] = linearizebmi(Fstr,...
                                                {Xstr,Ystr,'Z','M'},...
                                                {'X0dummy','Y0dummy','Z0dummy','M0dummy'},...
                                                '',lopts);
        catch
            warning("(linearizebmi) Execution failed, resize decision matricies and retry...")
            % resize Z and M
            [vstruct(vstep),Z,Z0dummy,M,M0dummy] = resizeZM(vstruct(vstep));
            % execute linearizebmi
            [LMIauto,~,gLMI,BMI] = linearizebmi(Fstr,...
                                                {Xstr,Ystr,'Z','M'},...
                                                {'X0dummy','Y0dummy','Z0dummy','M0dummy'},...
                                                '',lopts);
        end
    elseif isZ
        try
            [LMIauto,~,gLMI,BMI] = linearizebmi(Fstr,...
                                                {Xstr,Ystr,'Z'},...
                                                {'X0dummy','Y0dummy','Z0dummy'},...
                                                '',lopts);
        catch
            warning("(linearizebmi) Execution failed, resize decision matricies and retry...")
            % resize Z and M
            [vstruct(vstep),Z,Z0dummy] = resizeZM(vstruct(vstep));
            % execute linearizebmi
            [LMIauto,~,gLMI,BMI] = linearizebmi(Fstr,...
                                                {Xstr,Ystr,'Z'},...
                                                {'X0dummy','Y0dummy','Z0dummy'},...
                                                '',lopts);
        end
    else
        [LMIauto,~,gLMI,BMI] = linearizebmi(Fstr,...
                                            {Xstr,Ystr},...
                                            {'X0dummy','Y0dummy'},...
                                            '',lopts);
    end
    
    
    % LMIauto
    % gLMI.sdpvarname
    
    % declare LMI constraints
    if gLMI.isbmi
        if vlen > 1 && isempty(vlist{i})
            error("varargin{2}{%d} must not be empty cell list.\n\n(Corresponded BMI):\n%s",i,Fstr)
        end
        
        vstruct(vstep).dLMI = LMIauto;
        vstruct(vstep).BMIopt = gLMI;
        vstruct(vstep).orgBMI = BMI;
        if vstep < vlen
            vstep = vstep + 1;
        end
        bminum = bminum + 1;
    else
        if vlen > 1 && ~isempty(vlist{i}) 
            error("varargin{2}{%d} must be empty cell list '{}'.\n\n(Corresponded LMI):\n%s",i,Fstr)
        end
        LMIlist = [LMIlist, LMIauto<=-1e-6];
    end
    
    % sdpvar names in all constraints
    sdpvarnamelist = [sdpvarnamelist, gLMI.sdpvarname];
    
end

if ~(vstep == bminum)
    error("The number of BMIs (in varargin{1}) and set of decision variables (in varargin{2}) must be the same: \n(There are %d BMIs, %d vars set in input)",bminum,vlen);
end

% dLMI
% LMIlist
% sdpvarnamelist

% vstruct(1)
% vstruct(2)
% disp("line: 185")

sdpvarnamelist = unique(sdpvarnamelist,'stable');



%% search init value by optimize [LMI<t] (t->negative)

%%% epsilon
eps = 1e-6;
upeps = 1e3;
loweps = 1e-3;


%%% decision value 't'
% A dilated LMI constraint is below
% [Q  XN      [tI O
%  GY -G]  -   O  O]  <  O
%

t = sdpvar;   % objective function
LMIinit = []; % dilated LMI for searching initial solution
for i=1:bminum
    % size of linear term Q in the dilated LMI
    sizeQ = size(vstruct(i).BMIopt.data.Q);
    sizeLMI = size(vstruct(i).dLMI);

    % decide objective function 't'
    tI = t * eye(sizeQ);
    alpha = blkdiag(tI, zeros(sizeLMI-sizeQ));

    %%% constraint for search init val
    LMIi = [vstruct(i).dLMI <= alpha]; % minimize t
    if isZ
      LMIi = [LMIi, vstruct(i).Z+vstruct(i).Z'>=eps];
      %LMIi = [LMIi, Z+Z'>=loweps];
      %LMIi = [LMIi, Z+Z'<=upeps];
    end
    if isM
      LMIi = [LMIi, vstruct(i).M+vstruct(i).M'>=eps];    
    end
    
    LMIinit = [LMIinit, LMIi];
end

% Append other LMI
LMIinit = [LMIinit, LMIlist];

% (Just Test): Limits gamma bounds
LMIinit = [LMIinit, 1e3>=g>=0];


%%% Init val
% set init val by assign()
% if not assign, set default val
%
% All of these process is just for
% finding max eig of LMIauto as alpha(t)


% Decide others sdpvar default value
firstvaluelist = cell(1,length(sdpvarnamelist));
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
for i=1:bminum
    assign(vstruct(i).X0dummy,value(vstruct(i).X))
    assign(vstruct(i).Y0dummy,value(vstruct(i).Y))
    assign(vstruct(i).Z0dummy,eye(vstruct(i).sizeZ))
    assign(vstruct(i).M0dummy,eye(vstruct(i).sizeM))

    % set default value if not exist assign value
    X0init = value(vstruct(i).X0dummy);
    if isnan(X0init)
        % val = zeros(sizeX);
        val = eye(vstruct(i).sizeX);
        assign(vstruct(i).X0dummy,val)
        assign(vstruct(i).X,val)
        X0init = value(vstruct(i).X0dummy);
    end

    Y0init = value(vstruct(i).Y0dummy);
    if isnan(Y0init)
        val = zeros(vstruct(i).sizeY);
        assign(vstruct(i).Y0dummy,val)
        assign(vstruct(i).Y,val)
        Y0init = value(vstruct(i).Y0dummy);
    end

    Z0init = value(vstruct(i).Z0dummy);
    assign(vstruct(i).Z,Z0init)
    
    M0init = value(vstruct(i).M0dummy);
    assign(vstruct(i).M,M0init)
    
    vstruct(i).X0init = X0init;
    vstruct(i).Y0init = Y0init;
    vstruct(i).Z0init = Z0init;
    vstruct(i).M0init = M0init;
end


% Assign initial objective values (actually, not needed)
if isnan(value(g))
    assign(g,zeros(size(g)));
end


% Find max eig in dilated LMI with alpha
eiglmi = eig(value(blkdiag(vstruct(1:bminum).dLMI)));
maxeig = max(eiglmi);

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

maxeig

ttall = tt;   % optimal solutions
lc = 0;
while tt >= 0
  lc=lc+1;
  
  extLMI=LMIinit;
  for i=1:bminum
      extLMI=replace(extLMI,vstruct(i).X0dummy,vstruct(i).X0init);
      extLMI=replace(extLMI,vstruct(i).Y0dummy,vstruct(i).Y0init);
      extLMI=replace(extLMI,vstruct(i).Z0dummy,vstruct(i).Z0init);
      extLMI=replace(extLMI,vstruct(i).M0dummy,vstruct(i).M0init);
  end

  optimize(extLMI,t,opts.yalmip);
  for i=1:bminum
      vstruct(i).X0init=double(vstruct(i).X);
      vstruct(i).Y0init=double(vstruct(i).Y);
      vstruct(i).Z0init=double(vstruct(i).Z);
      vstruct(i).M0init=double(vstruct(i).M);
  end
  
  tt=double(t);
  ttall=[ttall,tt];
  
  % gout=double(g)
  
  % debug output (maybe not necessary)
  if debug
      fprintf('Loop#%03d: %9.4f\n',lc,tt)
  end
  
  
  % loop count upper bound
  if lc == 1000
      % disp("debug:loop:1")
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
LMI = [];
for i=1:bminum
    LMI = [LMI, vstruct(i).dLMI<=0];
    if isZ
      LMI = [LMI, vstruct(i).Z+vstruct(i).Z'>=eps];
      % LMI = [LMI, Z+Z'>=loweps];
      % LMI = [LMI, Z+Z'<=upeps];
    end
    if isZ
      LMI = [LMI, vstruct(i).M+vstruct(i).M'>=eps];
    end
end

% Append other LMI
LMI = [LMI, LMIlist];


% Add regularization term
v = sdpvar(1,1);
pt = opts.penalty;

if ~isequal(pt,0) && ~opts.regterm
    for i=1:bminum
        X = vstruct(i).X;
        Y = vstruct(i).Y;
        Z = vstruct(i).Z;
        M = vstruct(i).M;
        X0dummy = vstruct(i).X0dummy;
        Y0dummy = vstruct(i).Y0dummy;
        Z0dummy = vstruct(i).Z0dummy;
        M0dummy = vstruct(i).M0dummy;
        
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
        end
        
        if isM
          [vmr,vmc]=size(M);
          vm =sdpvar(vmr,vmr,'symmetric');
          if issymmetric(M)
            LMI=[LMI,[vm,triu(M-M0dummy);triu(M-M0dummy)',eye(vmc)]>=0];
          else
            LMI=[LMI,[vm,M-M0dummy;(M-M0dummy)',eye(vmc)]>=0];
          end
        end

        if isM
          LMI=[LMI,v>=trace(vx)+trace(vy)+trace(vz)+trace(vm)];
        elseif isZ
          LMI=[LMI,v>=trace(vx)+trace(vy)+trace(vz)];
        else
          LMI=[LMI,v>=trace(vx)+trace(vy)];
        end
    end
end


% init val from above process
for i=1:bminum
    vstruct(i).X0 = vstruct(i).X0init;
    vstruct(i).Y0 = vstruct(i).Y0init;
    if isZ
      vstruct(i).Z0 = vstruct(i).Z0init;
    end
    if isM
      vstruct(i).M0 = vstruct(i).M0init;
    end
end


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
  for i=1:bminum
      extLMI=replace(extLMI,vstruct(i).X0dummy,vstruct(i).X0);
      extLMI=replace(extLMI,vstruct(i).Y0dummy,vstruct(i).Y0);
      if isZ
        extLMI=replace(extLMI,vstruct(i).Z0dummy,vstruct(i).Z0);
      end
      if isM
        extLMI=replace(extLMI,vstruct(i).M0dummy,vstruct(i).M0);
      end
  end
  
  %%% applying regularization term
  if opts.regterm && opts.penalty>0
      lmdc = opts.penalty;
      terms = 0;
      for i=1:bminum
          terms = terms + regterm(lmdc,vstruct(i).X, vstruct(i).Y,...
                                       vstruct(i).X0,vstruct(i).Y0,lc);
      end
      optimize(extLMI,g+terms,opts.yalmip);
  else
      optimize(extLMI,g+v*pt,opts.yalmip);
  end
  
  %%% optimize by dilated LMI constraits as sufficient conditions of BMI
  % optimize(extLMI,g+v*pt,opts.yalmip);
  tEnd = toc(tStart);
  
  
  % update determined val
  for i=1:bminum
      vstruct(i).X0=double(vstruct(i).X);
      vstruct(i).Y0=double(vstruct(i).Y);
      if isZ
          vstruct(i).Z0=double(vstruct(i).Z);
      end
      if isM
          vstruct(i).M0=double(vstruct(i).M);
      end
  end

  
  
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
outopts.alphaall = ttall; % objective function for initial
% outopts.tmall = [0,tmall]; % computational time
outopts.tmall = tmall;


% optimal solution
for i=1:length(sdpvarnamelist)
    name = sdpvarnamelist(i);
    data = evalin('caller',name);
    vars.(name) = value(data);
end
%vars

% initial solution
for i=1:bminum
    vars.([vstruct(i).Xstr '_init']) = vstruct(i).X0init;
    vars.([vstruct(i).Ystr '_init']) = vstruct(i).Y0init;
    if isZ
      vars.(['G' num2str(i) '_init']) = vstruct(i).Z0init;
    end
    if isM
      vars.(['M' num2str(i) '_init']) = vstruct(i).M0init;
    end
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



%% update matrix Z and M when input {X,Y} are {Y,X} (inversion input order)
function [vstruct,Z,Z0dummy,M,M0dummy] = resizeZM(vstruct)

% size X
sizeX = vstruct.sizeX;
% update info about Z
if issymmetric(vstruct.Z)
    Z = sdpvar(sizeX(1),sizeX(1),'symmetric');
    sizeZ = size(Z);
    Z0dummy = sdpvar(sizeZ(1),sizeZ(2),'symmetric');
else
    Z = sdpvar(sizeX(1),sizeX(1),'full');
    sizeZ = size(Z);
    Z0dummy = sdpvar(sizeZ(1),sizeZ(2),'full');
end
vstruct.Z = Z;
vstruct.sizeZ = sizeZ;
vstruct.Z0dummy = Z0dummy;
% update info about M
M = sdpvar(sizeX(1),sizeX(1),'symmetric');
sizeM = size(M);
M0dummy = sdpvar(sizeM(1),sizeM(2),'symmetric');
vstruct.M = M;
vstruct.sizeM = sizeM;
vstruct.M0dummy = M0dummy;

