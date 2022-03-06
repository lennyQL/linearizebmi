%% H_2 static output feedback control problem
% To run this example, 
% required SDP solver 'SeDuMi'


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


%%% Define decision matrices
p=sdpvar(nx,nx,'symmetric');	% Lyapunov matrix
k=sdpvar(nu,ny,'full');         % Controller (static gain)
g=sdpvar(1,1);                  % H-inf norm
R=sdpvar(nz,nz,'symmetric');    % H2norm > sqrt(trace(R))

% Assign initial values to sdpvar
assign(p,eye(nx,nx))
assign(k,zeros(nu,ny))
assign(g,0)
assign(R,eye(nz,nz))


%%% YALMIP optimize settings
yalmipopts=sdpsettings;
yalmipopts.solver='sedumi';	% SDP solver
yalmipopts.verbose=0;       % No details

%%% solvebmi options:
opts = solvebmiOptions;
opts = solvebmiOptions(opts,'yalmip',yalmipopts);   % yalmip options
opts = solvebmiOptions(opts,'lcmax',200);            % Maximum step times
% opts.penalty=1e-2;

%%% BMI as a string
Fstr1 = "[p*(a+b2*k*c2)+(p*(a+b2*k*c2))',p*(b1+b2*k*d21);"+...
        "(p*(b1+b2*k*d21))',            -eye(nw,nw)];";

% others constraints
Fstr2 = "trace(R)-g";
Fstr3 = "[-R              -(c1+d12*k*c2);"+...
        "-(c1+d12*k*c2)'  -p];";
    
% All constraints (describe matrix variables)
% In solvebmi, "Fstr" can be define [Fstr<=0]
Flist = {Fstr1, Fstr2, Fstr3};


%%% Use solvebmi() to execute ovebounding aproximation method 
%%% to solve the BMI proble
%%% Use 2 dilation type to solve BMI

% (Type 1) Decomposition matrix G is a constant matrix
% :Solve speed is fast
opts = solvebmiOptions(opts,'method',0);
[gg,vars,output] = solvebmi(Flist,{'p','k'},g,opts);

% (Type 2) Decomposition matrix G is a decision matrix
% :Converges quickly
opts = solvebmiOptions(opts,'method',1);
[gg2,vars2,output2] = solvebmi(Flist,{'p','k'},g,opts);


%%% output
% Final acheived optimal value
gg
gg2

% Designed controller
K1 = vars.k
K2 = vars2.k

% Data for all defined sdpvar as optimal solutions
% vars
% vars2


%%% Figure achieved value g
ggall = [output.ggall; output2.ggall]';
figure;
plot(ggall,'LineWidth',2);
xlabel('Number of Iteration')
ylabel('$H_{\infty}$ norm','Interpreter','latex')
legend('dilated LMI (Sebe 2007)','dilated LMI (Sebe 2018)')
grid on


%%% Figure process about optimizing "alpha" 
%%% for searching initial feasible solutions
ttall = shapePlotData(output.ttall, output2.ttall);
figure;
plot(ttall,'LineWidth',2);
xlabel('Number of Iteration')
ylabel('$\alpha$','Interpreter','latex')
legend('dilated LMI (Sebe 2007)','dilated LMI (Sebe 2018)')
xticks(0:1:200)
grid on


