%%% SDP solver settings (YALMIP)
opts=sdpsettings;
opts.solver='sedumi';	% SDP solver
opts.verbose=0;
stoptol=5e-7;		% stop tolerance
lcmax=1e4;		% maximum step numbers
pt=1e-2;		% penalty factor in optimization


%%% Define problem (from COMPleib)
%[a,b1,b2,c1,c2,d11,d12,d21,nx,nw,nu,nz,ny]=COMPleib('AC3');
%------------------------------------------------------------------
%(AC3): L-1011 aircraft in cruise flight conditions
%       C. Edwards and S. K. Spurgeon,
%       On the development of discontinuous observers",
%       IJOC, Vol. 59, Nr. 5, pp. 1211-1229, 1994
%------------------------------------------------------------------
nx=5;nu=2;ny=4;
a=[0 0 1 0 0;0 -0.154 -0.0042 1.54 0;0 0.249 -1 -5.2 0;
   0.0386 -0.996 -0.0003 -0.117 0;0 0.5 0 0 -0.5];
b2=[0 0;-0.744 -0.032;0.337 -1.12;0.02 0;0 0];
c2=[0 1  0 0 -1;0 0 1 0 0;0 0 0 1 0;1 0 0 0 0];
b1=eye(nx); c1=eye(nx);
[nx,nw]=size(b1);
[nz,nx]=size(c1);
d12=[zeros(nz-nu,nu);eye(nu)];  
d11=zeros(nz,nw);
d21=zeros(ny,nw);

%% Definitions of decision matrices
% Define decision matrices
p=sdpvar(nx,nx,'symmetric');	% Lyapunov matrix
k=sdpvar(nu,ny,'full');         % Controller (static gain)
g=sdpvar(1,1);                  % H-infinity norm
v=sdpvar(1,1);			% penalty term

% Define dummy decision variables for successive overbounding
p0=sdpvar(nx,nx,'symmetric');
k0=sdpvar(nu,ny,'full');

% Define decomposition matrix
nk=size(k,1);
G =sdpvar(nk,nk,'full');   % decomposition matrix (decision value)
G0=sdpvar(nk,nk,'full');   % decomposition matrix (dummy)

%% Definition of BMI problem
% BMI as a string
Fstr = "[p*(a+b2*k*c2)+(p*(a+b2*k*c2))',p*(b1+b2*k*d21),(c1+d12*k*c2)';"   +...
        "(p*(b1+b2*k*d21))',            -g*eye(nw,nw),     (d11+d12*k*d21)';" +...
        "c1+d12*k*c2,                   d11+d12*k*d21,  -g*eye(nz)]";

    
%%% Use linearizebmi() to convert BMI into dilated LMI
% (case 1) decomposition matrix G is constant
%[LMIauto,LMIstr]=linearizebmi(Fstr,{'p','k'},{'p0','k0'},'G0');
% (case 2) decomposition matrix G is a decision value
[LMIauto,LMIstr,~,orgBMI]=linearizebmi(Fstr,{'p','k','G'},{'p0','k0','G0'});

% Define extended H-infinity constraints
LMI=[LMIauto<=0,p>=1e-6];

% add penalty term
vm=[vec(p)-vec(p0);vec(k)-vec(k0);vec(G)-vec(G0)];
vn=size(vm,1);
LMI=[LMI,[eye(vn),vm;vm',v]>=0];


%%% Overbounding Aprroximation Method
%% Initial feasible solutions
% (K=O is a stabilizing static gain)
k0init=zeros(nu,ny)	
%% Calculate initial Lyapunov matrix and H-infinity norm
initLMI=replace(orgBMI,k,k0init);
optimize([initLMI<=0,p>=1e-6],g,opts);
p0init=double(p)
ggsav=double(g)
G0init=eye(size(G));

%% sequentially optimization 
ggall=ggsav;
for lc=1:lcmax
  extLMI=LMI;
  extLMI=replace(extLMI,p0,p0init);
  extLMI=replace(extLMI,k0,k0init);
  extLMI=replace(extLMI,G0,G0init);

  optimize(extLMI,g+v*pt,opts);
%   optimize(extLMI,g,opts);

  p0init=double(p);
  k0init=double(k);
  G0init=double(G);

  % Print out achieved optimal value
  gg=double(g);
  ggall=[ggall,gg];
  fprintf('Loop#%03d: %9.4f\n',lc,gg);

  if ggsav-gg<stoptol
    break;
  end
  ggsav=gg;
end

%% Figures achieved H-infinity norm
% linear scale
figure;
plot(ggall);
grid on;

% log-scale
figure;
semilogy(ggall);
grid on;
