%%% SDP solver の設定
opts=sdpsettings;
opts.solver='sedumi';	% 使用する SDP solver
opts.verbose=0;

%%% 問題の定義
% 制御対象データの定義
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
% [a,b1,b2,c1,c2,d11,d12,d21,nx,nw,nu,nz,ny] = COMPleib('WEC1');
% [a,b1,b2,c1,c2,d11,d12,d21,nx,nw,nu,nz,ny] = COMPleib('WEC2');
% [a,b1,b2,c1,c2,d11,d12,d21,nx,nw,nu,nz,ny] = COMPleib('WEC3');


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
% 決定変数の定義
p=sdpvar(nx,nx,'symmetric');	% Lyapunov 行列
k=sdpvar(nu,ny,'full');		% 制御器(定数ゲイン)
g=sdpvar(1,1);			% H∞ノルム

% 双線形項決定変数の暫定解の宣言 (dummy)
% p0=sdpvar(nx,nx,'symmetric');
% k0=sdpvar(nu,ny,'full');
% 初期暫定解の宣言
p0=zeros(size(p));
k0=zeros(size(k));

% オプション:
% 最適化する値(need)
opts.g = g;
% 繰り返し実行する回数(need)
opts.lcmax = 200;
% 各ステップを表示(bool)
% opts.showstep = 0;
% 分割行列Gも最適化するか(bool)
% opts.dilate = 1;



% BMI 最適化問題の定義
Fstr = "[p*(a+b2*k*c2)+(p*(a+b2*k*c2))',p*(b1+b2*k*d21),(c1+d12*k*c2)';"   +...
        "(p*(b1+b2*k*d21))',            -g*eye(nw,nw),     (d11+d12*k*d21)';" +...
        "c1+d12*k*c2,                   d11+d12*k*d21,  -g*eye(nz)]";

%%% 提案した solvebmi() で 逐次 LMI 化法を実行
% 分割行列Gは定数
[gg,output] = solvebmi(Fstr,{'p','k'},opts);
% 分割行列Gは決定変数
opts.dilate = 1;
[gg2,output2] = solvebmi(Fstr,{'p','k'},opts);


gg
% output
% output.X
% output.Y


gg2
% output2
% output2.X
% output2.Y


eig(output.X)
eig(output2.X)

eig(output.X0)
eig(output2.X0)

eig(output2.Z+output2.Z')


% 達成値の更新過程の表示
% figure;
% plot(output.ggall);
% figure;
% semilogy(output.ggall);

% figure;
% plot(output2.ggall);
% figure;
% semilogy(output2.ggall);

ggall = [output.ggall; output2.ggall]';
figure;
plot(ggall);
figure;
semilogy(ggall);

