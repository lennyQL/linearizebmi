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
disp(newline)
disp("#####*** Hinf問題 ***#####")
% disp(newline)

%%% 決定変数の定義
p=sdpvar(nx,nx,'symmetric');	% Lyapunov 行列
k=sdpvar(nu,ny,'full');		% 制御器(定数ゲイン)
g=sdpvar(1,1);			% H∞ノルム
 

% 初期値
assign(p,eye(nx,nx))
assign(k,zeros(nu,ny))
assign(g,0)


%%% SDP solver の設定
yalmipopts=sdpsettings;
yalmipopts.solver='sedumi';	% 使用する SDP solver
yalmipopts.verbose=0;         % 詳細表示

%%% solvebmiのオプション:
opts.yalmip = yalmipopts;   % yalmipのoptimizeのためのオプション
opts.lcmax = 200;            % 繰り返し実行する回数
% opts.showstep = 0;        % 各ステップを表示(bool)
% opts.dilate = 1;          % 分割行列Gも最適化するか(bool)



%%% BMI 最適化問題の定義
Fstr = "[p*(a+b2*k*c2)+(p*(a+b2*k*c2))',p*(b1+b2*k*d21),(c1+d12*k*c2)';"   +...
        "(p*(b1+b2*k*d21))',            -g*eye(nw,nw),     (d11+d12*k*d21)';" +...
        "c1+d12*k*c2,                   d11+d12*k*d21,  -g*eye(nz)]";

    
% すべての制約
Flist = {Fstr, "-p"};


%%% 提案した solvebmi() で 逐次 LMI 化法を実行
%%%
% 分割行列Gは定数
opts.dilate = 0;
[gg,vars,output] = solvebmi(Flist,{'p','k'},g,opts);

%%%
% 分割行列Gは決定変数
opts.dilate = 1;
[gg2,vars2,output2] = solvebmi(Flist,{'p','k'},g,opts);


%%% output
% gg
% vars
% vars.g
% vars.p
% vars.k
% vars.p0
% vars.k0
% output.ttall'


% gg2
% vars2
% vars2.p
% vars2.k
% vars2.p0
% vars2.k0

% output2.ttall'

gg
gg2

K1 = vars.k
K2 = vars2.k

vars
vars2


%%% リアプノフ行列Pの正定値性の確認(固有値すべて正か)
% 最適解
% ispositive_g1_X = eig(output.X)
% ispositive_g2_X = eig(output2.X)

% 初期解
% ispositive_g1_X0 = eig(output.X0)
% ispositive_g2_X0 = eig(output2.X0)

% 分割行列
% ispositive_g2_Z = eig(output2.Z+output2.Z')


%%% 達成値の更新過程の表示
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
plot(ggall,'LineWidth',2);
xlabel('Number of Iteration')
ylabel('$H_{\infty}$ norm','Interpreter','latex')
legend('dilated LMI (5)','dilated LMI (7)')
grid on

% figure;
% semilogy(ggall);


%%% 時間経過
% ttall = [output.tmall; output2.tmall]';
% figure;
% plot(tmall,ggall,'LineWidth',2);
% xlabel('Computational Time')
% ylabel('$H_{\infty}$ norm','Interpreter','latex')
% legend('Sebe(2007)','Sebe(2018)')
% grid on


%%% alphaの最適化過程表示
[ttall1,ttall2] = matchSize(output.ttall, output2.ttall);
ttall = [ttall1; ttall2]';
figure;
plot(ttall,'LineWidth',2);
xlabel('Number of Iteration')
ylabel('$\alpha$','Interpreter','latex')
legend('dilated LMI (5)','dilated LMI (7)')
xticks(0:1:200)
grid on

% figure;
% semilogy(output2.ttall,output2.ggall);



%% H2問題
% solvebmiは現在，BMI制約を1つしか受け取れない
% その他のLMI制約をどう受け取るか
% 一応応急策として，optsに残りの制約を渡してる
disp(newline)
disp("#####*** H2問題 ***#####")
% disp(newline)

%%% 決定変数の定義
p=sdpvar(nx,nx,'symmetric');	% Lyapunov 行列
k=sdpvar(nu,ny,'full');         % 制御器(定数ゲイン)
g=sdpvar(1,1);                  % H2ノルム
R=sdpvar(nz,nz,'symmetric');    % H2ノルム=sqrt(trace(R))
 

% 初期値
assign(p,eye(nx,nx))
assign(k,zeros(nu,ny))
assign(R,zeros(nz,nz))
% assign(R,eye(nz,nz))


%%% SDP solver の設定
yalmipopts=sdpsettings;
yalmipopts.solver='sedumi';	% 使用する SDP solver
yalmipopts.verbose=0;         % 詳細表示
    
%%% solvebmiのオプション:
opts.yalmip = yalmipopts;   % yalmipのoptimizeのためのオプション
opts.lcmax = 200;            % 繰り返し実行する回数
% opts.showstep = 0;        % 各ステップを表示(bool)
% opts.dilate = 1;          % 分割行列Gも最適化するか(bool)



%%% BMI 最適化問題の定義
Fstr1 = "[p*(a+b2*k*c2)+(p*(a+b2*k*c2))',p*(b1+b2*k*d21);"+...
        "(p*(b1+b2*k*d21))',            -eye(nw,nw)];";

    
% すべての制約
Fstr2 = "trace(R)-g";
% 本来の[g*g>=trace(R)]は
% 非線形なのでsedumiで解けない
% テストとしてgの2乗部分を外した
% constraints = [g >= sqrtm(trace(R))];

Fstr3 = "[-R              -(c1+d12*k*c2);"+...
        "-(c1+d12*k*c2)'  -p];";

Flist = {Fstr1, Fstr2, Fstr3};


%%% 提案した solvebmi() で 逐次 LMI 化法を実行
%%%
% 分割行列Gは定数
opts.dilate = 0;
[gg,vars,output] = solvebmi(Flist,{'p','k'},g,opts);

%%%
% 分割行列Gは決定変数
opts.dilate = 1;
[gg2,vars2,output2] = solvebmi(Flist,{'p','k'},g,opts);


gg
gg2

K1 = vars.k
K2 = vars2.k

vars
vars2


%%% 達成値の更新過程の表示
ggall = [output.ggall; output2.ggall]';
figure;
plot(ggall,'LineWidth',2);
xlabel('Number of Iteration')
ylabel('$H_{2}$ norm','Interpreter','latex')
% legend('Sebe(2007)','Sebe(2018)')
legend('dilated LMI (5)','dilated LMI (7)')
grid on

%%% alphaの最適化過程表示
[ttall1,ttall2] = matchSize(output.ttall, output2.ttall);
ttall = [ttall1; ttall2]';
figure;
plot(ttall,'LineWidth',2);
xlabel('Number of Iteration')
ylabel('$\alpha$','Interpreter','latex')
% legend('Sebe(2007)','Sebe(2018)')
legend('dilated LMI (5)','dilated LMI (7)')
xticks(0:1:200)
grid on


%% 極配置問題(未完成)
% disp(newline)
% disp("#####*** 極配置問題 ***#####")
% disp(newline)

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



%% Useful Functions

function [M1,M2] = matchSize(L1,L2)
% 2つのリストのサイズを合わせる
% Nanを追加することで.
% グラフ出力用のデータ整形関数

sizeL1 = length(L1);
sizeL2 = length(L2);

M1 = L1;
M2 = L2;

sizediff = sizeL1 - sizeL2;

if sizediff > 0
    for i=1:abs(sizediff)
        M2 = [M2, NaN];
    end
elseif sizediff < 0
    for i=1:abs(sizediff)
        M1 = [M1, NaN];
    end
end


end


