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
%  T. Shimomura and T. Fujii,
%  Multiobjective control via successive over-bounding of quadratic terms.
%  International Journal of Robust and Nonlinear Control,
%  15, pp. 363–381, 2005.
%
% feedback system
a =[0,0,1,0;0,0,0,1;-5/4,5/4,0,0;5/4,-5/4,0,0];
b2=[0;0;1;0];
c2=[0,1,0,0];
d22=0;
nx=size(a,1);
nu=size(b2,2);
ny=size(c2,1);

% H2 performance
b0=[0;0;0;1];
c0=[0,1,0,0;0,0,0,0];
d00=[0;0];
d02=[0;2];
d20=1/10;
psys2=ss(a,[b0,b2],[c0;c2],[d00,d02;d20,d22]);

nw0=size(b0,2);
nz0=size(c0,1);

% H-infinity performance
b1=[0,0;0,0;-1/4,1/10;1/4,0];
c1=[1,-1,0,0;0,0,0,0];
d11=zeros(2,2);
d12=[0;1/5];
d21=[0,0];
psysinf=ss(a,[b1,b2],[c1;c2],[d11,d12;d21,d22]);
ginf=1;		% H-infinity norm bound

nw1=size(b1,2);
nz1=size(c1,1);


% order of controller
nk=nx;

% extended state space realization for nk~=0;
ea=blkdiag(zeros(nk),a);
eb2=blkdiag(eye(nk),b2);
ec2=blkdiag(eye(nk),c2);
ed22=blkdiag(zeros(nk),d22);

eb0=[zeros(nk,nw0);b0];
ec0=[zeros(nz0,nk),c0];
ed00=d00;
ed02=[zeros(nz0,nk),d02];
ed20=[zeros(nk,nw0);d20];

eb1=[zeros(nk,nw1);b1];
ec1=[zeros(nz1,nk),c1];
ed11=d11;
ed12=[zeros(nz1,nk),d12];
ed21=[zeros(nk,nw1);d21];



%%% SDP solver settings (YALMIP)
% initialization
yalmip('clear');

opts = solvebmiOptions('yalmip',sdpsettings);
opts.yalmip.verbose=0;
%
opts.yalmip.solver='sdpt3';	% SDP solver
%
% opts.yalmip.solver='sdpt3';	% SDP solver
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
%p2  =sdpvar(nx,nx,'symmetric');	% Lyapunov matrix (H2)
%pinf=sdpvar(nx,nx,'symmetric');	% Lyapunov matrix (Hinf)

% Lyapunov matrices (only for common Lyapunov approach)
x2  =sdpvar(nx,nx,'symmetric');	% Lyapunov matrix (H2)
y2  =sdpvar(nx,nx,'symmetric');
xinf=sdpvar(nx,nx,'symmetric');	% Lyapunov matrix (Hinf)
yinf=sdpvar(nx,nx,'symmetric');

% controller decision variables
ka=sdpvar(nx,nx,'full');
kb=sdpvar(nx,ny,'full');
kc=sdpvar(nu,nx,'full');
%kd=sdpvar(nu,ny,'full');
kd=zeros(nu,ny);
k=[ka,kb;kc,kd];
R=sdpvar(nz0,nz0,'symmetric');    % Complement matrix (H2)
g2=sdpvar(1,1);                 % H2 norm
ginf=1;                         % Hinf norm


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initial candidate of controller
% by common Lyapunov function approach.
% C. Scherer, P. Gahinet and M. Chilali,
% Multiobjective output-feedback control via LMI optimization.
% IEEE Transactions on Automatic Control, 42, pp. 896-911, 1997.

% common Lyapunov approach settings
LMI=[x2==xinf,y2==yinf];

% H2 constraints
M0=[a*x2+b2*kc,a,b0;...
    ka,y2*a+kb*c2,y2*b0+kb*d20;...
    zeros(nw0,nx+nx),-eye(nw0)/2];
M0=-(M0+M0');

M1=[x2,eye(nx),(c0*x2+d02*kc)';...
    eye(nx),y2,c0';...
    c0*x2+d02*kc,c0,R];

LMI=[LMI,M0>=0,M1>=0];

% H-infinity constraints
M2=[a*xinf+b2*kc,a+b2*kd*c2,b1+b2*kd*d21,zeros(nx,nz1);...
    ka,yinf*a+kb*c2,yinf*b1+kb*d21,zeros(nx,nz1);...
    zeros(nw1,nx+nx),-ginf/2*eye(nw1),zeros(nw1,nz1);...
    c1*xinf+d12*kc,c1+d12*kd*c2,d11+d12*kd*d21,-ginf/2*eye(nz1)];
M2=-(M2+M2');

LMI=[LMI,M2>=0];

optimize(LMI,trace(R),opts.yalmip);


% 
fprintf('Guaranteed H2 norm (common Lyapunov):  %9.6f\n',sqrt(trace(value(R))));

xv=value(x2);
yv=value(y2);
kav=value(ka);
kbv=value(kb);
kcv=value(kc);
kdv=kd;

[ku,ks,kv]=svd(eye(nx)-yv*xv);
ku=ku*sqrtm(ks);
kv=kv*sqrtm(ks);

p=[-kv\xv*ku,ku';ku,yv];
p=(p+p')/2;

kav=ku\(kav-[yv,kbv]*[a,b2;c2,zeros(ny,nu)]*[xv;kcv])/(kv');
kbv=ku\kbv;
kcv=kcv/(kv');

ksys=ss(kav,kbv,kcv,kdv);
kv=[kav,kbv;kcv,kdv];

%pole(lft(psys2,ksys))
fprintf('Achieved H2 norm (common Lyapunov):    %9.6f\n',norm(lft(psys2,  ksys),2,  1e-6));
fprintf('Achieved H-inf norm (common Lyapunov): %9.6f\n',norm(lft(psysinf,ksys),inf,1e-6));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% non-common Lyapunov function approach by solving BMI
% decision variables
p2  =sdpvar(nx*2,nx*2,'symmetric');	% Lyapunov matrix (H2)
pinf=sdpvar(nx*2,nx*2,'symmetric');	% Lyapunov matrix (Hinf)


% check code

M0=-[p2*(ea+eb2*k*ec2)+(p2*(ea+eb2*k*ec2))',p2*(eb0+eb2*k*ed20);...
    (p2*(eb0+eb2*k*ed20))',                -eye(nw0,nw0)];

M1=[R              (ec0+ed02*k*ec2);...
    (ec0+ed02*k*ec2)', p2];

M2=-[pinf*(ea+eb2*k*ec2)+(pinf*(ea+eb2*k*ec2))',pinf*(eb1+eb2*k*ed21),(ec1+ed12*k*ec2)';...
   (pinf*(eb1+eb2*k*ed21))',                -ginf*eye(nw1),    (ed11+ed12*k*ed21)';...
    ec1+ed12*k*ec2,                          ed11+ed12*k*ed21,    -ginf*eye(nz1)];

assign(p2,  p);
assign(pinf,p);
assign(k,kv);
%assign(R,value(R));

eig(value(M0))
eig(value(M1))
eig(value(M2))



% Describe BMIs and LMIs

Fstr1 = "[p2*(ea+eb2*k*ec2)+(p2*(ea+eb2*k*ec2))',p2*(eb0+eb2*k*ed20);"+...
        "(p2*(eb0+eb2*k*ed20))',                -eye(nw0)];";

Fstr2 = "[pinf*(ea+eb2*k*ec2)+(pinf*(ea+eb2*k*ec2))',pinf*(eb1+eb2*k*ed21),(ec1+ed12*k*ec2)';"+...
        "(pinf*(eb1+eb2*k*ed21))',                -ginf*eye(nw1),    (ed11+ed12*k*ed21)';"+...
        "ec1+ed12*k*ec2,                          ed11+ed12*k*ed21,    -ginf*eye(nz1)];";

%Fstr3 = "trace(R)-g2";

Fstr4 = "-[R              (ec0+ed02*k*ec2);"+...
        " (ec0+ed02*k*ec2)', p2];";

%Flist = {Fstr1, Fstr2, Fstr3, Fstr4,"-pinf"};
Flist = {Fstr1, Fstr2, Fstr4,"-pinf"};


    
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
  [gg{tc},vars{tc},output{tc}] = solvebmi(Flist,{{'p2','k'},{'pinf','k'}},trace(R),opts);
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
xticks(0:1:1000)
grid on;
hold off;

% Achieved controllers
for i=1:length(lgd)
  vars{i}.k
end
