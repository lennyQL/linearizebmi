%% 出力フィードバック制御のBMI制約についてのテスト
%   bmiparserを評価するテストプログラム
%   Hinf問題，極配置問題2つの制約それぞれ評価する

% 次元
n =4; % xの次元,状態数
m1=3; % wの次元
m2=2; % uの次元,入力数
p1=3; % zの次元
p2=2; % yの次元,出力数

% 決定変数の宣言
X=sdpvar(n,n, 'full');
Y=sdpvar(m2,p2, 'full');
% BMI 決定変数の暫定値
X0=rand(size(X));
Y0=rand(size(Y));

% 一般化制御対象の係数行列
% random例
A=rand(n,n);
B1=rand(n,m1);
B2=rand(n,m2);
C1=rand(p1,n); 
C2=rand(p2,n);
D11=rand(p1,m1);
D12=rand(p1,m2);
D21=rand(p2,m1);
D22=rand(p2,m2);

% G,gammaの定数倍
G = eye(size(Y,1));
% 単位行列
I = eye(p1);



%% Hinf問題, 制約がブロック行列の場合
disp("***** H無限大制御問題 *****")

%%%%%% parser(関数)による自動計算 %%%%%%%
% Hinf問題の制約の左辺
Fstr = "[X*A+X*B2*Y*C2+A'*X'+C2'*Y'*B2'*X'  X*B1+X*B2*Y*D21  C1'+C2'*Y'*D12';"+...
       "B1'*X'+D21'*Y'*B2'*X'               -I               D11'+D21'*Y'*D12';"+...
       "C1+D12*Y*C2                         D11+D12*Y*D21    -I]";
% parser(自作した関数)
[LMIauto, BMIauto] = bmiparser(Fstr,["X","Y"],["X0","Y0"]);


%%%%%% 手計算 %%%%%%
% BMI制約
BMImanual = [X*A+X*B2*Y*C2+A'*X'+C2'*Y'*B2'*X'  X*B1+X*B2*Y*D21  C1'+C2'*Y'*D12';
             B1'*X'+D21'*Y'*B2'*X'               -I               D11'+D21'*Y'*D12';
             C1+D12*Y*C2                         D11+D12*Y*D21    -I];
% 線形項,heあり
Q0 = [X*A+A'*X'    X*B1           C1'+C2'*Y'*D12';
      B1'*X'       -I             D11'+D21'*Y'*D12';
      C1+D12*Y*C2  D11+D12*Y*D21    -I];
% 双線形項の定数行列,heなし
L = [eye(n,n);zeros(m1,n);zeros(p1,n)];
N = B2;
R = [C2 D21 zeros(p2,p1)];
% 手計算で導いた逐次LMI,heあり
LMImanual = [Q0+L*X*N*Y0*R+L*X0*N*Y*R-L*X0*N*Y0*R+(L*X*N*Y0*R+L*X0*N*Y*R-L*X0*N*Y0*R)',...
             L*(X-X0)*N+R'*(Y-Y0)'*G';...
            (L*(X-X0)*N+R'*(Y-Y0)'*G')',...
             -2*G];


%%%%% 評価 %%%%%
BMItt = BMImanual - BMIauto
LMItt = LMImanual - LMIauto



%% 極配置問題, 制約がブロック行列で表されない場合
disp(newline)
disp("***** 極配置問題 *****")

% % BMI 制約式の左辺
Fstr = "X*A+X*B2*Y*C2+A'*X'+C2'*Y'*B2'*X'";
[LMIauto,BMIauto] = bmiparser(Fstr,["X","Y"],["X0","Y0"]);
% 
BMImanual = X*A+X*B2*Y*C2+A'*X'+C2'*Y'*B2'*X';
%
Q0 = X*A+A'*X';
L = eye(n,n);
N = B2;
R = C2;
%
LMImanual = [Q0+L*X*N*Y0*R+L*X0*N*Y*R-L*X0*N*Y0*R+(L*X*N*Y0*R+L*X0*N*Y*R-L*X0*N*Y0*R)',...
             L*(X-X0)*N+R'*(Y-Y0)'*G';...
            (L*(X-X0)*N+R'*(Y-Y0)'*G')',...
             -2*G];
%
BMItt = BMImanual - BMIauto
LMItt = LMImanual - LMIauto



