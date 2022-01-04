%%% initialization
yalmip('clear');

%%% SDP solver settings (YALMIP)
opts.yalmip=sdpsettings;
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

opts.max_retry=5;

% opts.method  
%   opts.method=0: Sebe (2007)
%   opts.method=1: Sebe (2018)
%   opts.method=2: Shimomura & Fujii (2005)
%   opts.method=3: Lee & Hu (2016)
%   opts.method=4: Ren et al. (2021)
opts.method1.t=0.99;
opts.method1.t=0;

stoptol=5e-7;		% stop tolerance
lcmax=2e2;		% maximum step numbers
regfac=1e-3;		% regularization factor in optimization


%%% Define problem (from COMPleib)
%probids={'AC3','AC6','AC15','AC16','AC17',...
%         'HE2','DIS1','DIS3','TG1','AGS','WEC2','WEC3',...
%         'MFP','UWV','PSM','NN4','NN8','NN11'};
%probid='NN4'; % 'AGS','MFP'
%
probid='NN4';
[a,b1,b2,c1,c2,d11,d12,d21,nx,nw,nu,nz,ny]=COMPleib(probid);

% DIS1: CS 系の方が良さそう
% DIS3: optimal への収束がばらつく
% TG1: CS 系の方が良さそう


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
G =sdpvar(nk,nk,'full');   % decision matrix
G0=sdpvar(nk,nk,'full');   % dummy matrix


% open loop H-infinity norm
F = [p*a+(p*a)',p*b1,c1'; ...
    (p*b1)',-g*eye(nw,nw),d11'; ...
     c1,d11,-g*eye(nz)]

optimize([F<=0,p>=1e-6],g,opts.yalmip);
p0init_g=double(p);
k0init_g=zeros(nu,ny);
double(g)

% BMI as a string
Fstr  = "[p*(a+b2*k*c2)+(p*(a+b2*k*c2))',p*(b1+b2*k*d21),(c1+d12*k*c2)';"   +...
        "(p*(b1+b2*k*d21))',            -g*eye(nw,nw),     (d11+d12*k*d21)';" +...
        "c1+d12*k*c2,                   d11+d12*k*d21,  -g*eye(nz)]";

    
%%% Use linearizebmi() to convert BMI into dilated LMI
%tall=[0.99,0.7,0.5,0,-1];
lgd={'Sebe (2007)','Sebe (2018) (t=0)','Sebe (2018) (t=0.99)',...
     'Shimoura & Fujii (2005)','Lee & Hu (2016)','Ren et al. (2021)'};
ggcomp=cell(length(lgd),1);
sgcomp=cell(length(lgd),1);
for tc=1:6
  % method
  if tc<3
    opts.method=tc-1;
  else
    opts.method=tc-2;
  end

  % for method 1
  if tc==3
    t=0.99;
  else
    t=0
  end
  if t<0
    F_t=1;
    t=sdpvar(1,1);
  else
    F_t=0;
  end

  opts.t = t;

  switch opts.method
    case 0  % Sebe (2007)
      G0=eye(nk);
      [LMIauto,LMIstr,~,orgBMI]...
          =linearizebmi(Fstr,{'p','k'},{'p0','k0'},'G0',opts);
    case 1  % Sebe (2018)
      G =sdpvar(nk,nk,'full');        % decomposition matrix
      G0=sdpvar(nk,nk,'full');        % dummy matrix
      [LMIauto,LMIstr,~,orgBMI]...
           =linearizebmi(Fstr,{'p','k','G'},{'p0','k0','G0'},'',opts);
    case 2  % Shimomura & Fujii (2005)
      G0=eye(nk);
      [LMIauto,LMIstr,~,orgBMI]...
          =linearizebmi(Fstr,{'p','k'},{'p0','k0'},'G0',opts);
    case 3  % Lee & Hu (2016)
      G =sdpvar(nk,nk,'symmetric');   % decomposition matrix
      G0=sdpvar(nk,nk,'symmetric');   % dummy matrix
      [LMIauto,LMIstr,~,orgBMI]...
           =linearizebmi(Fstr,{'p','k','G'},{'p0','k0','G0'},'',opts);
    case 4  % Ren el al. (2021)
      G =sdpvar(nk,nk,'full');        % decomposition matrix
      G0=sdpvar(nk,nk,'full');        % dummy matrix
      M =sdpvar(nk,nk,'symmetric');   % second decomposition matrix
      M0=sdpvar(nk,nk,'symmetric');   % dummy matrix
      [LMIauto,LMIstr,~,orgBMI]...
           =linearizebmi(Fstr,{'p','k','G','M'},{'p0','k0','G0','M0'},'',opts);
  end


  % Define constraints
  %  LMI=[LMIauto<=0,p>=1e-6];
  LMI=[LMIauto<=0,p>=1e-6];

  switch opts.method
    case {1,2,3,4}
      LMI=[LMI,[1e3*eye(nu),G;G',1e3*eye(nu)]>=0];
      LMI=[LMI,G+G'>=1e-3*eye(nu)];
    otherwise
      ;
  end
  switch opts.method
    case 4
      LMI=[LMI,[1e3*eye(nu),M;M',1e3*eye(nu)]>=0];
      LMI=[LMI,M+M'>=1e-3*eye(nu)];
    otherwise
      ;
  end
  

  % add reguralization term
  vm=[vec(p)-vec(p0);vec(k)-vec(k0);vec(G)-vec(G0)];
  vn=size(vm,1);
  LMI=[LMI,[eye(vn),vm;vm',v]>=0];


  %%% Overbounding Aprroximation Method
  % Initial feasible solutions
  % (K=O is a stabilizing static gain)
  k0init=zeros(nu,ny);	
  %% Calculate initial Lyapunov matrix and H-infinity norm
  initLMI=replace(orgBMI,k,k0init);
  optimize([initLMI<=0,p>=1e-6],g,opts.yalmip);
  p0init=double(p);
  ggsav=double(g)
  G0init=eye(size(G));
  M0init=eye(size(G));

  % sequentially optimization
  ggall=ggsav;
  sgall=max(svd(G0init));
  ggsav=ggsav*2;
  retry=0;
  for lc=1:lcmax
    extLMI=LMI;
    extLMI=replace(extLMI,p0,p0init);
    extLMI=replace(extLMI,k0,k0init);
    switch opts.method
      case {0,2}
        ;
      case {1,3}
        extLMI=replace(extLMI,G0,G0init);
      case 4
        extLMI=replace(extLMI,G0,G0init);
        extLMI=replace(extLMI,M0,M0init);
      otherwise
        ;
    end

    if F_t
      tlc=1-1/lc;
      tlc=min(0.99,max(1e-2,tlc));
      extLMI=replace(extLMI,t,tlc);
    end

    sol=optimize(extLMI,g+v*regfac,opts.yalmip);
    sol

    p0init=double(p);
    k0init=double(k);
    switch opts.method
      case {0,2}
        ;
      case {1,3}
        G0init=double(G)
      case 4
        G0init=double(G)
        M0init=double(M)
      otherwise
        ;
    end

    % Print out achieved optimal value
    gg=double(g);
    ggall=[ggall,gg];
    fprintf('Loop#%03d: %9.4f\n',lc,gg);

    sg=max(svd(G0init));
    sgall=[sgall,sg];

    %pause
    if ggsav-gg<stoptol
      retry=retry+1;
      if retry==opts.max_retry
        break;
      end
    else
      retry=0;
      ggsav=gg;
      ksav =k0init;
    end

  end
  ggcomp{tc}=ggall;
  sgcomp{tc}=sgall;
  kcomp{tc} =ksav;

end

figure;
for i=1:length(ggcomp)
  semilogy(0:length(ggcomp{i})-1,ggcomp{i});
  hold on;
end
grid on;
legend(lgd);
title(probid);

figure;
for i=1:length(sgcomp)
  semilogy(0:length(sgcomp{i})-1,sgcomp{i});
  hold on;
end
grid on;
legend(lgd);
title(probid);

for i=1:length(kcomp)
  kcomp{i}
end
