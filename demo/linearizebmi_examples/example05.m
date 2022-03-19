%%% ---------------------------------------------------------- %%%
% H-infinity control by static output feedback
% with the overbounding approximation proposed in Sebe (2018).
% 
% Additionally,
%  - append regularization terms to objective function.
%  - compare with several parameters 't' for decomposition matrix.
%%% ---------------------------------------------------------- %%%



%%% SDP solver settings (YALMIP)
opts=sdpsettings;
opts.verbose=0;
%
opts.solver='sedumi';	% SDP solver
%
% opts.solver='sdpt3';	% SDP solver
opts=sdpsettings(opts,'sdpt3.gaptol',1e-8);
opts=sdpsettings(opts,'sdpt3.inftol',1e-10);
opts=sdpsettings(opts,'sdpt3.steptol',1e-10);
opts=sdpsettings(opts,'sdpt3.maxit',5e2);
opts=sdpsettings(opts,'sdpt3.predcorr',1);
opts=sdpsettings(opts,'sdpt3.expon',1);
opts=sdpsettings(opts,'sdpt3.scale_data',0);
opts=sdpsettings(opts,'sdpt3.stoplevel',0);	% changed from default
opts=sdpsettings(opts,'sdpt3.printyes',0);

stoptol=5e-7;		% stop tolerance
lcmax=2e2;		% maximum step numbers
regfac=1e-4;		% regularization factor in optimization


%%% Define problem (from COMPleib)
%probids={'AC3','AC6','AC15','AC16','AC17',...
%         'AGS','DIS1','DIS3','HE2','MFP',...
%         'NN4','NN8','NN11','PSM','TG1',...
%         'UWV','WEC2','WEC3'};
%
probid='MFP';
[a,b1,b2,c1,c2,d11,d12,d21,nx,nw,nu,nz,ny]=COMPleib(probid);


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

optimize([F<=0,p>=1e-6],g,opts);
p0init_g=double(p);
k0init_g=zeros(nu,ny);
double(g)

% BMI as a string
Fstr  = "[p*(a+b2*k*c2)+(p*(a+b2*k*c2))',p*(b1+b2*k*d21),(c1+d12*k*c2)';"   +...
        "(p*(b1+b2*k*d21))',            -g*eye(nw,nw),     (d11+d12*k*d21)';" +...
        "c1+d12*k*c2,                   d11+d12*k*d21,  -g*eye(nz)]";

    
%%% Use linearizebmi() to convert BMI into dilated LMI
tall=[0.99,0.7,0.5,0,-1];
lgd={'t=0.99','t=0.7','t=0.5','t=0','t=1/lc'};
ggcomp=cell(length(tall),1);
sgcomp=cell(length(tall),1);
for tc=1:length(tall)
  t=tall(tc);
  if t<0
    F_t=1;
    t=sdpvar(1,1);
  else
    F_t=0;
  end

  options = linearizebmiOptions('t',t,'method',1);
  % or:
  % options = linearizebmiOptions;
  % options.t = t;
  % options.method = 1;
  [LMIauto,LMIstr,~,orgBMI]...
      =linearizebmi(Fstr,{'p','k','G'},{'p0','k0','G0'},'',options);

  % Define constraints
  %  LMI=[LMIauto<=0,p>=1e-6];
  LMI=[LMIauto<=0,p>=1e-6,[1e5*eye(nu),G;G',1e5*eye(nu)]>=0];
  vm=[vec(p)-vec(p0);vec(k)-vec(k0);vec(G)-vec(G0)];
  vn=size(vm,1);
  LMI=[LMI,[eye(vn),vm;vm',v]>=0];


  %%% Overbounding Aprroximation Method
  % Initial feasible solutions
  % (K=O is a stabilizing static gain)
  k0init=zeros(nu,ny);	
  %% Calculate initial Lyapunov matrix and H-infinity norm
  initLMI=replace(orgBMI,k,k0init);
  optimize([initLMI<=0,p>=1e-6],g,opts);
  p0init=double(p);
  ggsav=double(g)
  G0init=eye(size(G));

  % sequentially optimization
  ggall=ggsav;
  sgall=max(svd(G0init));
  ggsav=ggsav*2;
  for lc=1:lcmax
    extLMI=LMI;
    extLMI=replace(extLMI,p0,p0init);
    extLMI=replace(extLMI,k0,k0init);
    extLMI=replace(extLMI,G0,G0init);
    if F_t
      tlc=1-1/lc;
      tlc=min(0.99,max(1e-2,tlc));
      extLMI=replace(extLMI,t,tlc);
    end

    sol=optimize(extLMI,g+v*regfac,opts);
    % sol

    p0init=double(p);
    k0init=double(k);
    G0init=double(G);

    % Print out achieved optimal value
    gg=double(g);
    ggall=[ggall,gg];
    fprintf('Loop#%03d: %9.4f\n',lc,gg);

    sg=max(svd(G0init));
    sgall=[sgall,sg];

    %pause
    if ggsav-gg<stoptol
      break;
    end
    ggsav=gg;

  end
  ggcomp{tc}=ggall;
  sgcomp{tc}=sgall;

end

figure;
for i=1:length(tall);
  semilogy(0:length(ggcomp{i})-1,ggcomp{i});
  hold on;
end
grid on;
legend(lgd);

figure;
for i=1:length(tall);
  semilogy(0:length(sgcomp{i})-1,sgcomp{i});
  hold on;
end
grid on;
legend(lgd);
