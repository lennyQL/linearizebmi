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
% [a,b1,b2,c1,c2,d11,d12,d21,nx,nw,nu,nz,ny] = COMPleib('AC17');
% [a,b1,b2,c1,c2,d11,d12,d21,nx,nw,nu,nz,ny] = COMPleib('WEC1');
% [a,b1,b2,c1,c2,d11,d12,d21,nx,nw,nu,nz,ny] = COMPleib('WEC2');
% [a,b1,b2,c1,c2,d11,d12,d21,nx,nw,nu,nz,ny] = COMPleib('WEC3');
% % [a,b1,b2,c1,c2,d11,d12,d21,nx,nw,nu,nz,ny] = COMPleib('MFP');


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


% nx=5;nu=2;ny=4;
% a=[0 0 1 0 0;0 -0.154 -0.0042 1.54 0;0 0.249 -1 -5.2 0;
%    0.0386 -0.996 -0.0003 -0.117 0;0 0.5 0 0 -0.5];
% b2=[0 0;-0.744 -0.032;0.337 -1.12;0.02 0;0 0];
% c2=[0 1  0 0 -1;0 0 1 0 0;0 0 0 1 0;1 0 0 0 0];
% b1=eye(nx); c1=eye(nx);
% [nx,nw]=size(b1);
% [nz,nx]=size(c1);
% d12=[zeros(nz-nu,nu);eye(nu)];  
% d11=zeros(nz,nw);
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

%%% デフォルトオプション
opts = solvebmiOptions;

%%% SDP solver の設定
yalmipopts=sdpsettings;
yalmipopts.solver='sdpt3';	% 使用する SDP solver
yalmipopts.verbose=0;         % 詳細表示

%%% solvebmiのオプション:
opts.yalmip = yalmipopts;   % yalmipのoptimizeのためのオプション
opts.lcmax = 200;            % 繰り返し実行する回数



%%% BMI 最適化問題の定義
Fstr = "[p*(a+b2*k*c2)+(p*(a+b2*k*c2))',p*(b1+b2*k*d21),(c1+d12*k*c2)';"   +...
        "(p*(b1+b2*k*d21))',            -g*eye(nw,nw),     (d11+d12*k*d21)';" +...
        "c1+d12*k*c2,                   d11+d12*k*d21,  -g*eye(nz)]";

    
% すべての制約
Flist = {Fstr, "-p"};


%%% 提案した solvebmi() で 逐次 LMI 化法を実行
%%% 分割行列のタイプ別の比較
% % 分割行列Gは定数
% opts.method = 0;
% opts.regterm= 0;
% [gg,vars,output] = solvebmi(Flist,{'p','k'},g,opts);
% %%%
% % 分割行列Gは決定変数
% opts.method = 1;
% opts.regterm= 0;
% [gg2,vars2,output2] = solvebmi(Flist,{'p','k'},g,opts);

%%% ペナルティ項ありなしの比較(分割行列は決定変数)
% ペナルティ項なし
opts.method = 1;
opts.penalty= 0;
opts.test = 0;
opts.testg = 0;
% [gg,vars,output] = solvebmi(Flist,{'p','k'},g,opts);
% ペナルティ項あり
opts.method = 1;
opts.penalty= 0;
opts.test = 0;
opts.testg = 1;
% [gg2,vars2,output2] = solvebmi(Flist,{'p','k'},g,opts);
% ペナルティ項あり
opts.method = 1;
opts.penalty= 0;
opts.test = 1;
opts.testg = 0;
% [gg3,vars3,output3] = solvebmi(Flist,{'p','k'},g,opts);
% ペナルティ項あり
opts.method = 1;
opts.penalty= 1e-2;
opts.test = 0;
opts.testg = 0;
% [gg4,vars4,output4] = solvebmi(Flist,{'p','k'},g,opts);
% ペナルティ項あり
opts.method = 1;
opts.penalty= 1e-2;
opts.test = 0;
opts.testg = 1;
% [gg5,vars5,output5] = solvebmi(Flist,{'p','k'},g,opts);
% ペナルティ項あり
opts.method = 1;
opts.penalty= 1e-2;
opts.test = 1;
opts.testg = 0;
% [gg6,vars6,output6] = solvebmi(Flist,{'p','k'},g,opts);

%%%
opts.method = 0;
opts.penalty= 0;
opts.test = 0;
opts.testg = 1;
[gg2,vars2,output2] = solvebmi(Flist,{'p','k'},g,opts);
%
opts.method = 1;
opts.penalty= 0;
opts.test = 0;
opts.testg = 1;
[gg3,vars3,output3] = solvebmi(Flist,{'p','k'},g,opts);
%
opts.method = 1;
opts.penalty= 1e-2;
opts.test = 0;
opts.testg = 0;
% [gg5,vars5,output5] = solvebmi(Flist,{'p','k'},g,opts);
%
opts.method = 1;
opts.penalty= 1e-4;
opts.test = 0;
opts.testg = 1;
% [gg6,vars6,output6] = solvebmi(Flist,{'p','k'},g,opts);



%%% output
% gg
% vars
% vars.g
% vars.p
% vars.k
% vars.p0
% vars.k0
% output.alphaall'


% gg2
% vars2
% vars2.p
% vars2.k
% vars2.p0
% vars2.k0

% output2.alphaall'

% %
% gg
% gg2
% 
% K1 = vars.k
% K2 = vars2.k
% 
% vars
% vars2


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
% ggall = shapePlotData(output.ggall,output2.ggall,output3.ggall,output4.ggall,output5.ggall,output6.ggall)
% ggall = shapePlotData(output2.ggall,output3.ggall,output5.ggall,output6.ggall)
% ggall = shapePlotData(output2.ggall,output3.ggall,output6.ggall)
ggall = shapePlotData(output2.ggall,output3.ggall)
figure;
% plot(ggall,'LineWidth',1.5);
semilogy(ggall,'LineWidth',1.5);
xlabel('Number of Iteration')
ylabel('$H_{\infty}$ norm','Interpreter','latex')
% legend('Sebe(2007)','Sebe(2018)')
legend('Sebe(2007)','Sebe(2018)','Sebe(2018) : penalty')
% legend('None','Penalty')
% ylim([4 1e2])
grid on

% figure;
% semilogy(ggall);


%%% 時間経過
% tmall = shapePlotData(output2.tmall,output3.tmall)
% % tmall = shapePlotData(output2.tmall,output3.tmall,output5.tmall,output6.tmall)
% figure;
% % plot(tmall,ggall,'LineWidth',1.5);
% semilogy(tmall,ggall,'LineWidth',1.5);
% xlabel('Computational Time')
% ylabel('$H_{\infty}$ norm','Interpreter','latex')
% legend('Sebe(2007)','Sebe(2018)')
% % legend('None    : Alpha','None    : K0','Penalty : Alpha','Penalty : K0')
% grid on


%%% alphaの最適化過程表示
alphaall = shapePlotData(output2.alphaall,output3.alphaall)
% alphaall = shapePlotData(output2.alphaall,output3.alphaall,output6.alphaall)
figure;
plot(alphaall,'LineWidth',1.5);
xlabel('Number of Iteration')
ylabel('$\alpha$','Interpreter','latex')
legend('Sebe(2007)','Sebe(2018)','Sebe(2018) : penalty')
% legend('Sebe(2007)','Sebe(2018)')
xticks(0:1:200)
grid on





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
% opts.method = 1;          % 分割行列Gも最適化するか(bool)



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
%%% 分割行列のタイプ別の比較
% % 分割行列Gは定数
% opts.method = 0;
% opts.regterm= 0;
% [gg,vars,output] = solvebmi(Flist,{'p','k'},g,opts);
% %%%
% % 分割行列Gは決定変数
% opts.method = 1;
% opts.regterm= 0;
% [gg2,vars2,output2] = solvebmi(Flist,{'p','k'},g,opts);

%%% ペナルティ項ありなしの比較(分割行列は決定変数)
% ペナルティ項なし
opts.method = 1;
opts.regterm= 0;
[gg,vars,output] = solvebmi(Flist,{'p','k'},g,opts);
% ペナルティ項あり
opts.method = 1;
opts.regterm= 1;
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
plot(ggall,'LineWidth',1);
xlabel('Number of Iteration')
ylabel('$H_{2}$ norm','Interpreter','latex')
% legend('Sebe(2007)','Sebe(2018)')
legend('methodd LMI (5)','methodd LMI (7)')
grid on

%%% alphaの最適化過程表示
[alphaall1,alphaall2] = matchSize(output.alphaall, output2.alphaall);
alphaall = [alphaall1; alphaall2]';
figure;
plot(alphaall,'LineWidth',1);
xlabel('Number of Iteration')
ylabel('$\alpha$','Interpreter','latex')
% legend('Sebe(2007)','Sebe(2018)')
legend('methodd LMI (5)','methodd LMI (7)')
xticks(0:1:200)
grid on


%% 極配置問題(調整中)
disp(newline)
disp("#####*** 極配置問題 ***#####")
disp(newline)

% 決定変数の定義
p=sdpvar(nx,nx,'symmetric');	% Lyapunov 行列
k=sdpvar(nu,ny,'full');         % 制御器(定数ゲイン)
g=sdpvar(1,1);                  % 最大極の上限

% kとgを1つの決定変数に
b = [b2 -eye(size(b2,1))];
c = [c2; eye(size(c2,2))];
ka = blkdiag(k,g*eye(size(a)));


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


% 逐次LMI
Fstr = "p*(a+b*ka*c)+(p*(a+b*ka*c))'";
Flist = {Fstr,"-p"};

%%% 提案した solvebmi() で 逐次 LMI 化法を実行
%%% 分割行列のタイプ別の比較
% % 分割行列Gは定数
% opts.method = 0;
% opts.regterm= 0;
% [gg,vars,output] = solvebmi(Flist,{'p','ka'},g,opts);
% %%%
% % 分割行列Gは決定変数
% opts.method = 1;
% opts.regterm= 0;
% [gg2,vars2,output2] = solvebmi(Flist,{'p','ka'},g,opts);

%%% ペナルティ項ありなしの比較(分割行列は決定変数)
% ペナルティ項なし
opts.method = 1;
opts.regterm= 0;
[gg,vars,output] = solvebmi(Flist,{'p','ka'},g,opts);
% ペナルティ項あり
opts.method = 1;
opts.regterm= 1;
[gg2,vars2,output2] = solvebmi(Flist,{'p','ka'},g,opts);


% 出力
gg
gg2

K1 = vars.ka(1:nu,1:ny)
K2 = vars2.ka(1:nu,1:ny)

vars
vars2



%%% 達成値の更新過程の表示
ggall = [output.ggall; output2.ggall]';
figure;
plot(ggall,'LineWidth',1);
xlabel('Number of Iteration')
ylabel('maximum eigenvalue','Interpreter','latex')
legend('methodd LMI (5)','methodd LMI (7)')
grid on

%%% alphaの最適化過程表示
[alphaall1,alphaall2] = matchSize(output.alphaall, output2.alphaall);
alphaall = [alphaall1; alphaall2]';
figure;
plot(alphaall,'LineWidth',1);
xlabel('Number of Iteration')
ylabel('$\alpha$','Interpreter','latex')
legend('methodd LMI (5)','methodd LMI (7)')
xticks(0:1:200)
grid on



%% H2/Hinf問題
disp(newline)
disp("#####*** H2/Hinf問題 ***#####")


%%% 決定変数の定義
pinf=sdpvar(nx,nx,'symmetric');	% Lyapunov 行列 (Hinf)
p2=sdpvar(nx,nx,'symmetric');	% Lyapunov 行列 (H2)
k=sdpvar(nu,ny,'full');         % 制御器(定数ゲイン)
g=sdpvar(1,1);                  % H2ノルム
R=sdpvar(nz,nz,'symmetric');    % H2ノルム=sqrt(trace(R))
 

% 初期値
assign(pinf,eye(nx,nx))
assign(p2,eye(nx,nx))
assign(k,zeros(nu,ny))
% assign(R,zeros(nz,nz))
assign(R,eye(nz,nz))
% assign(g,0)


%%% SDP solver の設定
yalmipopts=sdpsettings;
yalmipopts.solver='sdpt3';	% 使用する SDP solver
yalmipopts.verbose=0;         % 詳細表示
    
%%% solvebmiのオプション:
opts = solvebmiOptions('yalmip',yalmipopts,... % yalmipのoptimizeのためのオプション
                       'lcmax', 200);          % 繰り返し実行する回数    


%%% BMI 最適化問題の定義
Fstr1 = "[p2*(a+b2*k*c2)+(p2*(a+b2*k*c2))',p2*(b1+b2*k*d21);"+...
        "(p2*(b1+b2*k*d21))',              -eye(nw,nw)];";

Fstr2 = "[pinf*(a+b2*k*c2)+(pinf*(a+b2*k*c2))',pinf*(b1+b2*k*d21),(c1+d12*k*c2)';"   +...
        "(pinf*(b1+b2*k*d21))',                -eye(nw,nw),       (d11+d12*k*d21)';" +...
        "c1+d12*k*c2,                          d11+d12*k*d21,     -eye(nz)]";

% その他のLMI
Fstr3 = "trace(R)-g";

Fstr4 = "[-R              -(c1+d12*k*c2);"+...
        "-(c1+d12*k*c2)'  -p2];";

Flist = {Fstr1, Fstr2, Fstr3, Fstr4,"-p2","-pinf"};


%%% 関数仕様test
opts.method = 0;
opts.penalty= 0;
% solvebmi(Flist,{{'p2','k'},{'pinf','k'},{'pinf','k'}},g,opts);
% [gg,vars,output] = solvebmi(Flist,{'pinf','k'},g,opts);
% [gg,vars,output] = solvebmi(Flist,{{'p2','k'},{'pinf','k'}},g,opts);
% ペナルティ項あり
opts.method = 4;
opts.penalty= 1e-6;
[gg2,vars2,output2] = solvebmi(Flist,{{'p2','k'},{'pinf','k'}},g,opts);


gg
gg2

K1 = vars.k
K2 = vars2.k

vars
vars2


% %%% 達成値の更新過程の表示
ggall = shapePlotData(output.ggall,output2.ggall)
figure;
% plot(ggall,'LineWidth',1);
semilogy(ggall,'LineWidth',1);
xlabel('Number of Iteration')
ylabel('$H_{\infty}$ norm','Interpreter','latex')
legend('Sebe(2007)','Sebe(2018)')
grid on

% %%% alphaの最適化過程表示
alphaall = shapePlotData(output.alphaall,output2.alphaall)
figure;
plot(alphaall,'LineWidth',1);
xlabel('Number of Iteration')
ylabel('$\alpha$','Interpreter','latex')
legend('Sebe(2007)','Sebe(2018)')
xticks(0:1:200)
grid on



