%% 極配置で初期解を決める

%%% SDP solver の設定
opts=sdpsettings;
opts.solver='sedumi';	% 使用する SDP solver
opts.verbose=0;

%%% 制御対象の係数行列，次元をインポート
[a,b1,b2,c1,c2,d11,d12,d21,nx,nw,nu,nz,ny] = COMPleib('HE1'); % g減るのに極増 
% [a,b1,b2,c1,c2,d11,d12,d21,nx,nw,nu,nz,ny] = COMPleib('HE2'); % 最初から極負
% [a,b1,b2,c1,c2,d11,d12,d21,nx,nw,nu,nz,ny] = COMPleib('HE3'); % 大,遅,HE4と似てる?
% [a,b1,b2,c1,c2,d11,d12,d21,nx,nw,nu,nz,ny] = COMPleib('HE4'); % lc40でギリギリ極負
% [a,b1,b2,c1,c2,d11,d12,d21,nx,nw,nu,nz,ny] = COMPleib('NN4'); % 最初から極負
% [a,b1,b2,c1,c2,d11,d12,d21,nx,nw,nu,nz,ny] = COMPleib('NN9'); % ランク欠落?
% [a,b1,b2,c1,c2,d11,d12,d21,nx,nw,nu,nz,ny] = COMPleib('NN17'); % size小,g振動
% [a,b1,b2,c1,c2,d11,d12,d21,nx,nw,nu,nz,ny] = COMPleib('AC3'); % 最初から極負
% [a,b1,b2,c1,c2,d11,d12,d21,nx,nw,nu,nz,ny] = COMPleib('AC7'); % lc:4で極負
%%% 初期解を変えたら一応解決

% a=[-0.0366,  0.0271,  0.0188, -0.4555;...
%     0.0482, -1.01,    0.0024, -4.0208;...
%     0.1002,  0.3681, -0.7070,  1.42;...
%     0,       0,       1,       0];
% b1=[4.678e-2, 0; 4.572e-2, 9.88e-3; 4.369e-2, 1.11e-3; -2.179e-2, 0];
% b2=[0.4422 0.1761;3.5446 -7.5922;-5.52 4.49;0 0];
% c1=(1/sqrt(2))*[2,0,0,0;0,1,0,0];
% c2=[0 1 0 0];
% nx=size(a,1);
% nw=size(b1,2);
% nu=size(b2,2);
% nz=size(c1,1);
% ny=size(c2,1);
% d11=zeros(nz,nw);
% d12=(1/sqrt(2))*[1, 0; 0, 1];
% d21=zeros(ny,nw);

%%

% 可制御性チェック
Co = ctrb(a,b2)              % フルランクか
unco = length(a) - rank(Co)  % 非可制御状態数

% 安定性チェック（極）
St = eig(a)


% 決定変数の定義
p=sdpvar(nx,nx,'symmetric');	% Lyapunov 行列
k=sdpvar(nu,ny,'full');         % 制御器(定数ゲイン)
z=sdpvar(nu,nu);      	% 分割行列
g=sdpvar(1,1);                  % H∞ノルム

% dummy
p0=sdpvar(nx,nx,'symmetric');	% Lyapunov 行列
k0=sdpvar(nu,ny,'full');		% 制御器(定数ゲイン)
z0=sdpvar(nu,nu);
t = sdpvar;


%%
%%% test
% g = sdpvar;
% Fstr = "[p*(a+b2*k*c2)+(p*(a+b2*k*c2))',p*(b1+b2*k*d21),(c1+d12*k*c2)';"   +...
%         "(p*(b1+b2*k*d21))',            -g*eye(nw,nw),     (d11+d12*k*d21)';" +...
%         "c1+d12*k*c2,                   d11+d12*k*d21,  -g*eye(nz)]";
% [LMIauto,LMIstr] = linearizebmi(Fstr,{'p','k'},{'p0','k0'})    
%%%%


% 極配置問題のための式変形, k,gを1つの決定変数とする（構造持ち:対角要素にk,gを並べる）
b = [b2 -eye(size(b2,1))];
c = [c2; eye(size(c2,2))];
ka = blkdiag(k,g*eye(size(a)));
ka0 = sdpvar(size(ka,1),size(ka,2));
% ka0 = blkdiag(sdpvar(nu,ny,'full'),-sdpvar*eye(size(a)));

% value(ka)


% kaのsizeに合わせて分割行列も拡大する，構造をどうするかの疑問は残る 
za = sdpvar(size(ka,1),size(ka,1));
% za = blkdiag(z,zeros(size(a)));
% za = blkdiag(z,eye(size(a)));
% za=sdpvar(nu,nu,'full');
za0 =sdpvar(size(za,1));




% BMIの線形化
Feig = "p*(a+b*ka*c)+(p*(a+b*ka*c))'";

% g = p;
% ka=k,ka0=k0,za=z,za0=z;
% Feig = "p*(a+b2*ka*c2)+(p*(a+b2*ka*c2))'-2*g";

% [LMIauto,LMIstr] = linearizebmi(Feig,{'p','ka'},{'p0','ka0'});        % Zなし
[LMIauto,LMIstr] = linearizebmi(Feig,{'p','ka','za'},{'p0','ka0','za0'}) % Zあり


% 許容範囲
eps = 1e-6;
% 制約の許容範囲をどうするか
% LMI = [LMIauto<=eps*eye(size(LMIauto))];
LMI = [LMIauto<=-eps*eye(size(LMIauto))];
% LMI = [LMIauto<=0];
% LMI = [LMI, p>=eps*eye(size(p))];


%%% run: overbounding approximation method
p0init = eye(size(p));
% k0の初期値をI*max(pole)
k0init = 1.1*max(real(eig(a))) * blkdiag(zeros(size(k)),eye(size(p)))
% z0init = zeros(size(za));
z0init = eye(size(za));

%

lcmax=200;
% ggall=[1.1*max(real(eig(a)))];
% maxeigall =[max(real(eig(a)))]; 
ggall=[];
maxeigall=[];

%%% 本来は極がすべて負になった時点で停止してよい
% maxeig = 1e3;
% lc = 0;
% while maxeig >= 0
%   lc=lc+1;
for lc=1:lcmax
  extLMI=LMI;
  extLMI=replace(extLMI,p0,p0init);
  extLMI=replace(extLMI,ka0,k0init);
  extLMI=replace(extLMI,za0,z0init);

  optimize(extLMI,g,opts);

  p0init=double(p);
  k0init=double(ka);
  z0init=double(za);
  
  % 制御器を抽出
  k0 = k0init(1:size(k,1),1:size(k,2));

  % show each step optimized value
  fbsys = a+b2*k0*c2;
  fbeig = real(eig(fbsys))
%   max(fbeig)
%   min(fbeig)
  maxeig=max(fbeig);
  maxeigall=[maxeigall,maxeig];

  gg=double(g);
  ggall=[ggall,gg];
  
  if ~isfield(opts,'showstep') || opts.showstep
      fprintf('Loop#%03d: %9.4f  %9.4f\n',lc,gg,maxeig);
%       fprintf('Loop#%03d: %9.4f\n',lc,gg)
%       fprintf('Loop#%03d: \n',lc);
%       disp(gg);
  end
  
  
  % ループ数上限到達で停止
  if lc == 30
      break
  end
  
end

fbsys = a+b2*k0*c2;
fbeig = real(eig(fbsys))
max(fbeig)
min(fbeig)

k0

lc

p0init
k0init
z0init

allplot = [ggall; maxeigall]';
% allplot = [maxeigall];

% 達成値の更新過程の表示
figure;
plot(allplot);
% figure;
% semilogy(allplot);

%%
% オプション:
% 最適化する値(need)
opts.g = g;
% 繰り返し実行する回数(need)
opts.lcmax = 200;
% 各ステップを表示
% opts.showstep = false;

% 極配置で求めた初期解
p0 = p0init;
k0 = k0init(1:size(k,1),1:size(k,2));
% z0 = z0init(1:size(z,1),1:size(z,2));
% p0=zeros(size(p));
% k0=zeros(size(k));
z0= eye(size(z));


% BMI 最適化問題の定義
Fstr = "[p*(a+b2*k*c2)+(p*(a+b2*k*c2))',p*(b1+b2*k*d21),(c1+d12*k*c2)';"   +...
        "(p*(b1+b2*k*d21))',            -g*eye(nw,nw),     (d11+d12*k*d21)';" +...
        "c1+d12*k*c2,                   d11+d12*k*d21,  -g*eye(nz)]";

%%% 提案した solvebmi() で 逐次 LMI 化法を実行
[gg,output] = solvebmi(Fstr,{'p','k'},{'p0','k0'},opts);

%%% 分割行列も決定変数の場合
[gg2,output2] = solvebmi(Fstr,{'p','k','z'},{'p0','k0','z0'},opts);


gg
output
output.X
output.Y

gg2
output2
output2.X
output2.Y
output2.Z


% 達成値の更新過程の表示
% figure;
% plot(output.ggall);
% figure;
% semilogy(output.ggall);

ggall = [output.ggall; output2.ggall]';
figure;
plot(ggall);
figure;
semilogy(ggall);


%% symsで式変換テスト

% syms a b2 c2 p k g I
% 
% b = [b2 1]
% c = [c2;1]
% ka = [k 0; 0 -g]
% 
% f = p*(a+b*ka*c)
% simplify(f)







