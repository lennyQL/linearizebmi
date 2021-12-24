%%% SDP solver settings (YALMIP)
opts=sdpsettings;
opts.solver='sedumi';	% SDP solver
opts.verbose=0;

%%% Define problem
% Define plant data
%  REFFERENCE:
%  S. N. Singh and A. A. R. Coelho,
%  "Nonlinear control of mismatched uncertain linear systems
%   and application to control of aircraft",
%  Journal of Dynamic Systems, Measurement and Control, 106,
%  pp.203-210, 1984
a=[-0.0366,  0.0271,  0.0188, -0.4555;...
    0.0482, -1.01,    0.0024, -4.0208;...
    0.1002,  0.3681, -0.7070,  1.42;...
    0,       0,       1,       0];
b1=[4.678e-2, 0; 4.572e-2, 9.88e-3; 4.369e-2, 1.11e-3; -2.179e-2, 0];
b2=[0.4422 0.1761;3.5446 -7.5922;-5.52 4.49;0 0];
c1=(1/sqrt(2))*[2,0,0,0;0,1,0,0];
c2=[0 1 0 0];
nx=size(a,1);
nw=size(b1,2);
nu=size(b2,2);
nz=size(c1,1);
ny=size(c2,1);
d11=zeros(nz,nw);
d12=(1/sqrt(2))*[1, 0; 0, 1];
d21=zeros(ny,nw);


% Define decision matrices
p=sdpvar(nx,nx,'symmetric');	% Lyapunov matrix
k=sdpvar(nu,ny,'full');         % Controller (static gain)
g=sdpvar(1,1);                  % H-inf norm


% Define feasible solutions for 
% each decision matrix in a bilinear term.
p0=sdpvar(nx,nx,'symmetric');
k0=sdpvar(nu,ny,'full');


% Define decomposition matrix
G =sdpvar(size(k,1));   % decision matrix
G0=eye(size(G));        % feasible solution


% BMI as a string
Fstr = "[p*(a+b2*k*c2)+(p*(a+b2*k*c2))',p*(b1+b2*k*d21),(c1+d12*k*c2)';"   +...
        "(p*(b1+b2*k*d21))',            -g*eye(nw,nw),     (d11+d12*k*d21)';" +...
        "c1+d12*k*c2,                   d11+d12*k*d21,  -g*eye(nz)]";

    
%%% Use linearizebmi() to convert BMI into dilated LMI
% (case 1) G is a constant matrix
[LMIauto,LMIstr]=linearizebmi(Fstr,{'p','k'},{'p0','k0'});
% (case 2) G is a decision matrix
% [LMIauto,LMIstr]=linearizebmi(Fstr,{'p','k','G'},{'p0','k0','G0'});

% Define constraints
LMI=[LMIauto<=0];


%%% Overbounding Aprroximation Method
% Initial feasible solutions
p0init=zeros(size(p));
k0init=zeros(size(k));
G0init=eye(size(G));

% sequentially optimization
lcmax=200;	% maximum step numbers
ggall=[];
for lc=1:lcmax
  extLMI=LMI;
  extLMI=replace(extLMI,p0,p0init);
  extLMI=replace(extLMI,k0,k0init);
  extLMI=replace(extLMI,k0,k0init);

  optimize(extLMI,g,opts);

  p0init=double(p);
  k0init=double(k);

  % Print out achieved optimal value
  gg=double(g);
  ggall=[ggall,gg];
  fprintf('Loop#%03d: %9.4f\n',lc,gg);

end

% Figure achieved value g
figure;
plot(ggall);
figure;
semilogy(ggall);

