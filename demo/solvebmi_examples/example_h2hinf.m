%%% ---------------------------------------------------------- %%%
% H_2/H_inf static output feedback control problem
% To run this example, 
% required SDP solver 'SDPT3'
% 
% Used method (research):
%     1. Sebe (2007)
%        - decomposition matrix 'G' is a constant matrix
%        
%     2. Sebe (2018)
%        - decomposition matrix 'G' is a decision matrix
%        - change parameter 't' in decomposition matrix
%     
%     3. Shimomura & Fujii (2005)
%        - overbounding approximation by 
%          completing the square with constant matrix
%     
%     4. Lee & Hu (2016)
%        - overbounding approximation by 
%          completing the square with decision matrix
%     
%     5. Ren et al. (2021)
%        - combining Sebe (2007) and Shimomura & Fujii (2005) 
%          for decomposition matrix
%%% ---------------------------------------------------------- %%%



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

% Or, use data set from COMPleib
probid='HE1';
% [a,b1,b2,c1,c2,d11,d12,d21,nx,nw,nu,nz,ny]=COMPleib(probid);


%%% initialization
yalmip('clear');


%%% SDP solver settings (YALMIP)
opts = solvebmiOptions('yalmip',sdpsettings);
opts.yalmip.verbose=0;
%
% opts.yalmip.solver='sedumi';	% SDP solver
%
opts.yalmip.solver='sdpt3';	% SDP solver
opts.yalmip=sdpsettings(opts.yalmip,'sdpt3.gaptol',1e-8);
opts.yalmip=sdpsettings(opts.yalmip,'sdpt3.inftol',1e-10);
opts.yalmip=sdpsettings(opts.yalmip,'sdpt3.steptol',1e-10);
opts.yalmip=sdpsettings(opts.yalmip,'sdpt3.maxit',5e2);
opts.yalmip=sdpsettings(opts.yalmip,'sdpt3.predcorr',1);
opts.yalmip=sdpsettings(opts.yalmip,'sdpt3.expon',1);
opts.yalmip=sdpsettings(opts.yalmip,'sdpt3.scale_data',0);
opts.yalmip=sdpsettings(opts.yalmip,'sdpt3.stoplevel',0);	% changed from default
opts.yalmip=sdpsettings(opts.yalmip,'sdpt3.printyes',0);

%%% Settings (solvebmi)
% opts.method  
%   opts.method=0: Sebe (2007)
%   opts.method=1: Sebe (2018)
%   opts.method=2: Shimomura & Fujii (2005)
%   opts.method=3: Lee & Hu (2016)
%   opts.method=4: Ren et al. (2021)
%
opts = solvebmiOptions(opts,'stoptol',5e-7);    % stop tolerance
opts = solvebmiOptions(opts,'lcmax',2e2);		% maximum step numbers
opts = solvebmiOptions(opts,'penalty',1e-4);	% regularization factor in optimization



%% Definitions of decision matrices
% Define decision matrices
p2=sdpvar(nx,nx,'symmetric');	% Lyapunov matrix (H2)
pinf=sdpvar(nx,nx,'symmetric');	% Lyapunov matrix (Hinf)
k=sdpvar(nu,ny,'full');         % Controller
R=sdpvar(nz,nz,'symmetric');    % Complement matrix (H2)
g2=sdpvar(1,1);                 % H2 norm
ginf=1;                         % Hinf norm

% Assign initial solution
assign(pinf,eye(nx,nx))
assign(p2,eye(nx,nx))
assign(k,zeros(nu,ny))
assign(R,eye(nz,nz))
assign(g2,0)



% Describe BMIs and LMIs
Fstr1 = "[p2*(a+b2*k*c2)+(p2*(a+b2*k*c2))',p2*(b1+b2*k*d21);"+...
        "(p2*(b1+b2*k*d21))',              -eye(nw,nw)];";

Fstr2 = "[pinf*(a+b2*k*c2)+(pinf*(a+b2*k*c2))',pinf*(b1+b2*k*d21),(c1+d12*k*c2)';"   +...
        "(pinf*(b1+b2*k*d21))',                -ginf*eye(nw,nw),  (d11+d12*k*d21)';" +...
        "c1+d12*k*c2,                          d11+d12*k*d21,     -ginf*eye(nz)]";

Fstr3 = "trace(R)-g2";

Fstr4 = "[-R              -(c1+d12*k*c2);"+...
        "-(c1+d12*k*c2)'  -p2];";

Flist = {Fstr1, Fstr2, Fstr3, Fstr4,"-p2","-pinf"};


    
%%% Use solvebmi() to execute ovebounding aproximation method 
lgd={'Sebe (2007)',...
     'Sebe (2018)',...
     'Shimoura & Fujii (2005)',...
     'Lee & Hu (2016)',...
     'Ren et al. (2021)'};

gg=cell(length(lgd),1);
vars=cell(length(lgd),1);
output=cell(length(lgd),1);
for tc=1:length(lgd)
  fprintf('\n--------------------\n');
  fprintf('Method: %s\n\n',lgd{tc});
  
  % Set method
  opts = solvebmiOptions(opts,'method',tc-1);
  
  % Execute solvebmi:
  % Decision variables       {{'p2','k'},{'pinf','k'}}
  % must be corresponding to {Fstr1,Fstr2}
  [gg{tc},vars{tc},output{tc}] = solvebmi(Flist,{{'p2','k'},{'pinf','k'},{},{},{},{}},g2,opts);
end


%%
% figure: H2 norm
figure;
for i=1:length(lgd)
  semilogy(0:length(output{i}.ggall)-1,output{i}.ggall);
  hold on;
end
legend(lgd);
xlabel('Number of Iteration')
ylabel('$H_{2}$ norm','Interpreter','latex')
title(probid);
grid on;
hold off;

% figure: alpha for searching initial solution
figure;
for i=1:length(lgd)
  plot(0:length(output{i}.alphaall)-1,output{i}.alphaall);
  hold on;
end
legend(lgd);
xlabel('Number of Iteration')
ylabel('$\alpha$','Interpreter','latex')
title(probid);
grid on;
hold off;

% Achieved controllers
for i=1:length(lgd)
  vars{i}.k
end
