%%% 問題の定義
% 制御対象データの定義
%%% 制御対象の係数行列，次元をインポート
[a,b1,b2,c1,c2,d11,d12,d21,nx,nw,nu,nz,ny] = COMPleib('HE1'); 
% [a,b1,b2,c1,c2,d11,d12,d21,nx,nw,nu,nz,ny] = COMPleib('HE2');
% [a,b1,b2,c1,c2,d11,d12,d21,nx,nw,nu,nz,ny] = COMPleib('HE3');
% [a,b1,b2,c1,c2,d11,d12,d21,nx,nw,nu,nz,ny] = COMPleib('HE4'); 
% [a,b1,b2,c1,c2,d11,d12,d21,nx,nw,nu,nz,ny] = COMPleib('NN4');
% [a,b1,b2,c1,c2,d11,d12,d21,nx,nw,nu,nz,ny] = COMPleib('NN9'); 
% [a,b1,b2,c1,c2,d11,d12,d21,nx,nw,nu,nz,ny] = COMPleib('NN17'); 
% [a,b1,b2,c1,c2,d11,d12,d21,nx,nw,nu,nz,ny] = COMPleib('AC3'); 
% [a,b1,b2,c1,c2,d11,d12,d21,nx,nw,nu,nz,ny] = COMPleib('AC7'); 
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


%% Hinf問題

%%% 決定変数の定義
p=sdpvar(nx,nx,'symmetric');	% Lyapunov 行列
k=sdpvar(nu,ny,'full');		% 制御器(定数ゲイン)
g=sdpvar(1,1);			% H∞ノルム



%%% SDP solver の設定
opts=sdpsettings;
opts.solver='sedumi';	% 使用する SDP solver
opts.verbose=0;         % 詳細表示

%%% solvebmiのオプション:
% opts.lcmax = 200;     % 繰り返し実行する回数
% opts.showstep = 0;    % 各ステップを表示(bool)
% opts.dilate = 1;      % 分割行列Gも最適化するか(bool)



%%% BMI 最適化問題の定義
Fstr = "[p*(a+b2*k*c2)+(p*(a+b2*k*c2))',p*(b1+b2*k*d21),(c1+d12*k*c2)';"   +...
        "(p*(b1+b2*k*d21))',            -g*eye(nw,nw),     (d11+d12*k*d21)';" +...
        "c1+d12*k*c2,                   d11+d12*k*d21,  -g*eye(nz)]";

%%% 提案した solvebmi() で 逐次 LMI 化法を実行
% 分割行列Gは定数
[gg,output] = solvebmi(Fstr,{'p','k'},g);
% 分割行列Gは決定変数
opts.dilate = 1;
[gg2,output2] = solvebmi(Fstr,{'p','k'},g,opts);


gg
% output
% output.X
% output.Y


gg2
% output2
% output2.X
% output2.Y


%%% リアプノフ行列Pの正定値性の確認(固有値すべて正か)
% 最適解
ispositive_g1_X = eig(output.X)
ispositive_g2_X = eig(output2.X)

% 初期解
ispositive_g1_X0 = eig(output.X0)
ispositive_g2_X0 = eig(output2.X0)

% 分割行列
ispositive_g2_Z = eig(output2.Z+output2.Z')


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



%% 極配置問題

% % 決定変数の定義
% p=sdpvar(nx,nx,'symmetric');	% Lyapunov 行列
% k=sdpvar(nu,ny,'full');         % 制御器(定数ゲイン)
% g=sdpvar(1,1);                  % 最大極の上限
% 
% % kとgを1つの決定変数に
% b = [b2 -eye(size(b2,1))];
% c = [c2; eye(size(c2,2))];
% ka = blkdiag(k,g*eye(size(a)));
% 
% 
% % オプション:
% % 最適化する値(need)
% % opts.g = g;
% % 繰り返し実行する回数(need)
% % opts.lcmax = 200;
% % 各ステップを表示(bool)
% % opts.showstep = 0;
% 
% 
% % 逐次LMI
% Fstr = "p*(a+b*ka*c)+(p*(a+b*ka*c))'";
% % 分割行列Gは定数
% opts.dilate = 0;
% [gg,output] = solvebmi(Fstr,{'p','ka'},g,opts);
% % 分割行列Gは決定変数
% opts.dilate = 1;
% [gg2,output2] = solvebmi(Fstr,{'p','ka'},g,opts);
% 
% 
% % 出力
% gg
% 
% gg2
% 
% output.Y
% output2.Y
% 
% eig(output.X)
% eig(output2.X)
% 
% eig(output.X0)
% eig(output2.X0)
% 
% eig(output2.Z+output2.Z')
% 
% 
% % 達成値の更新過程の表示
% ggall = [output.ggall; output2.ggall]';
% figure;
% plot(ggall);
% figure;
% semilogy(ggall);


