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
% ゼロ行列
O = zeros(p1);


% テスト用テーブルの初期化
cll = {};


%% Hinf問題, 制約がブロック行列の場合
disp("<<<======***** H無限大制御問題 *****======>>>")

%%%%%% 手動計算，真値 %%%%%%
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
             -(G+G')];
         

         
%%%%%% parser(関数)による自動計算 %%%%%%%
% すべて等価な表現
%% パターン1: 括弧()なし，1つの行列
disp(newline)
disp("###*** パターン1 ***###")
disp("# 括弧()なし，1つの行列による記述")

% Hinf問題の制約の左辺，bmiparserのための制約の記述
Fstr = "[X*A+X*B2*Y*C2+A'*X'+C2'*Y'*B2'*X'  X*B1+X*B2*Y*D21  C1'+C2'*Y'*D12';"+...
        "B1'*X'+D21'*Y'*B2'*X'               -I               D11'+D21'*Y'*D12';"+...
        "C1+D12*Y*C2                         D11+D12*Y*D21    -I]";

% parser(自作した関数)
tStart = tic;
[LMIauto, BMIauto] = bmiparser(Fstr,{'X','Y'},{'X0','Y0'});
tEnd = toc(tStart)

% 評価，真値との差
% BMItt: ゼロ行列かどうか
% LMItt: coefficient range < 1e-10 かどうか or ゼロ行列かどうか
BMItt = BMImanual - BMIauto;
LMItt = LMImanual - LMIauto;

disp("BMItt: "+evaltest(BMItt));
disp("LMItt: "+evaltest(LMItt));

cll = cat(1,cll,{"1. 括弧なし",evaltest(BMItt),evaltest(LMItt),tEnd});
    

    
%% パターン2: 括弧()あり，1つの行列
disp(newline)
disp("###*** パターン2 ***###")
disp("# 括弧()あり，1つの行列による記述")

Fstr = "[X*(A+B2*Y*C2)+(A+B2*Y*C2)'*X'  X*(B1+B2*Y*D21)  (C1+D12*Y*C2)';"+...
        "(B1+B2*Y*D21)'*X'               -I               (D11+D12*Y*D21)';"+...
        "C1+D12*Y*C2                     D11+D12*Y*D21    -I]";
    
tStart = tic;
[LMIauto, BMIauto] = bmiparser(Fstr,{'X','Y'},{'X0','Y0'});
tEnd = toc(tStart)

BMItt = BMImanual - BMIauto;
LMItt = LMImanual - LMIauto;

disp("BMItt: "+evaltest(BMItt));
disp("LMItt: "+evaltest(LMItt));

cll = cat(1,cll,{"2. 括弧あり",evaltest(BMItt), evaltest(LMItt),tEnd});
    
    
%% パターン3: 関数あり，1つの行列
disp(newline)
disp("###*** パターン3 ***###")
disp("# 関数あり，1つの行列")

Fstr = "[X*(A+B2*Y*C2)+(A+B2*Y*C2)'*X'  X*(B1+B2*Y*D21)     (C1+D12*Y*C2)';"+...
        "(B1+B2*Y*D21)'*X'               -eye(p1)+zeros(p1)  (D11+D12*Y*D21)';"+...
        "C1+D12*Y*C2                     D11+D12*Y*D21       -eye(p1,p1)*eye(p1)]";
    
tStart = tic;
[LMIauto, BMIauto] = bmiparser(Fstr,{'X','Y'},{'X0','Y0'});
tEnd = toc(tStart)

BMItt = BMImanual - BMIauto;
LMItt = LMImanual - LMIauto;

disp("BMItt: "+evaltest(BMItt));
disp("LMItt: "+evaltest(LMItt));

cll = cat(1,cll,{"3. 関数あり",evaltest(BMItt),evaltest(LMItt),tEnd});
    
    
%% パターン4: 複数の行列の和
disp(newline)
disp("###*** パターン4 ***###")
disp("# 複数の行列の和")

Fstr = "[[X*(A+B2*Y*C2)  X*(B1+B2*Y*D21)  zeros(n,p1);"+...
        "zeros(m1,n)     zeros(p1,m1)     zeros(m1,p1);"+...
        "zeros(p1,n)     zeros(p1,m1)     zeros(p1,p1)]"...     % 双線形項
        + "+" +...
        "[(A+B2*Y*C2)'*X'   zeros(n,m1)   zeros(n,p1);"+...
        "(B1+B2*Y*D21)'*X'  zeros(p1,m1)  zeros(m1,p1);"+...
        "zeros(p1,n)        zeros(p1,m1)  zeros(p1,p1)]"...     % 双線形項の転置
        + "+" +...
        "[zeros(n,n)   zeros(n,m1)    (C1+D12*Y*C2)';"+...
        "zeros(m1,n)   -eye(p1)       (D11+D12*Y*D21)';"+...
        "C1+D12*Y*C2   D11+D12*Y*D21  -eye(p1)]]";               % 線形項
    
tStart = tic;
[LMIauto, BMIauto] = bmiparser(Fstr,{'X','Y'},{'X0','Y0'});
tEnd = toc(tStart)

BMItt = BMImanual - BMIauto;
LMItt = LMImanual - LMIauto;

disp("BMItt: "+evaltest(BMItt));
disp("LMItt: "+evaltest(LMItt));

cll = cat(1,cll,{"4. 行列同士の和",evaltest(BMItt),evaltest(LMItt),tEnd});

    
%% パターン5: ベクトル同士の積, 転置あり
disp(newline)
disp("###*** パターン5 ***###")
disp("# ベクトル同士の積")

Fstr = "[[X;zeros(m1,n);zeros(p1,n)]*[A+B2*Y*C2 B1+B2*Y*D21 zeros(n,p1)]"... % 双線形項
        +"+"+...
        "([X;zeros(m1,n);zeros(p1,n)]*[A+B2*Y*C2 B1+B2*Y*D21 zeros(n,p1)])'"... % 双線形項の転置
        + "+" +...
        "[zeros(n,n)   zeros(n,m1)    (C1+D12*Y*C2)';"+...
        "zeros(m1,n)   -eye(p1)       (D11+D12*Y*D21)';"+...
        "C1+D12*Y*C2   D11+D12*Y*D21  -eye(p1)]]";            % 線形項
    
tStart = tic;
[LMIauto, BMIauto] = bmiparser(Fstr,{'X','Y'},{'X0','Y0'});
tEnd = toc(tStart)

BMItt = BMImanual - BMIauto;
LMItt = LMImanual - LMIauto;

disp("BMItt: "+evaltest(BMItt));
disp("LMItt: "+evaltest(LMItt));

cll = cat(1,cll,{"5. ベクトル同士の積",evaltest(BMItt),evaltest(LMItt),tEnd});
    
    
%% パターン6
disp(newline)
disp("###*** パターン6 ***###")
disp("# 構造の異なる入れ子表現1")

Fstr = "[[X;zeros(m1,n);zeros(p1,n)]*[A+B2*Y*C2 B1+B2*Y*D21 zeros(n,p1)]"... % 双線形項
        +"+"+...
        "([X;zeros(m1,n);zeros(p1,n)]*[A+B2*Y*C2 B1+B2*Y*D21 zeros(n,p1)])'"... % 双線形項の転置
        + "+" +...
        "[zeros(n,p1);zeros(m1,p1);eye(p1)]*[C1+D12*Y*C2 D11+D12*Y*D21 zeros(p1)]"... % 線形項
        + "+" +...
        "([zeros(n,p1);zeros(m1,p1);eye(p1)]*[C1+D12*Y*C2 D11+D12*Y*D21 zeros(p1)])'"... % 線形項の転置
        + "+" +...
        "-blkdiag(zeros(n),eye(m1),eye(p1))]";  % 対角ブロックの単位行列
    
tStart = tic;
[LMIauto, BMIauto] = bmiparser(Fstr,{'X','Y'},{'X0','Y0'});
tEnd = toc(tStart)

BMItt = BMImanual - BMIauto;
LMItt = LMImanual - LMIauto;

disp("BMItt: "+evaltest(BMItt));
disp("LMItt: "+evaltest(LMItt));

cll = cat(1,cll,{"6. 入れ子表現1",evaltest(BMItt),evaltest(LMItt),tEnd});

    
%% パターン7
disp(newline)
disp("###*** パターン7 ***###")
disp("# 構造の異なる入れ子表現2")

Fstr = "[[X*(A+B2*Y*C2)+(A+B2*Y*C2)'*X'  X*(B1+B2*Y*D21)     (C1+D12*Y*C2)';"+...
        "(B1+B2*Y*D21)'*X'               zeros(p1)           (D11+D12*Y*D21)';"+...
        "C1+D12*Y*C2                     D11+D12*Y*D21       zeros(p1)]"...
        + "+" +...
        "blkdiag(zeros(n),-eye(m1+p1))+zeros(10)]";
    
tStart = tic;
[LMIauto, BMIauto] = bmiparser(Fstr,{'X','Y'},{'X0','Y0'});
tEnd = toc(tStart)

BMItt = BMImanual - BMIauto;
LMItt = LMImanual - LMIauto;

disp("BMItt: "+evaltest(BMItt));
disp("LMItt: "+evaltest(LMItt));

cll = cat(1,cll,{"7. 入れ子表現2",evaltest(BMItt),evaltest(LMItt),tEnd});


%% パターン8
disp(newline)
disp("###*** パターン8 ***###")
disp("# blkbmiを使った表現")

% blkbmiのインスタント
Fstr = blkbmi; 
% ブロック行列の各要素での宣言
% 対称部分(転置)は記述しなくてよい
Fstr(1,1) = "X*(A+B2*Y*C2)+(A+B2*Y*C2)'*X'";
Fstr(1,2) = "X*(B1+B2*Y*D21)";
Fstr(2,2) = "-eye(p1)";
Fstr(3,1) = "C1+D12*Y*C2";
Fstr(3,2) = "D11+D12*Y*D21";
Fstr(3,3) = "-eye(p1,p1)";
    
tStart = tic;
[LMIauto, BMIauto] = bmiparser(Fstr,{'X','Y'},{'X0','Y0'});
tEnd = toc(tStart)

BMItt = BMImanual - BMIauto;
LMItt = LMImanual - LMIauto;

disp("BMItt: "+evaltest(BMItt));
disp("LMItt: "+evaltest(LMItt));

cll = cat(1,cll,{"8. blkbmiによる表現",evaltest(BMItt),evaltest(LMItt),tEnd});


%% 極配置問題, 制約がブロック行列で表されない場合
disp(newline)
disp("<<<======***** 極配置問題 *****======>>>")


% 手動計算，真値
BMImanual = X*(A+B2*Y*C2)+(A+B2*Y*C2)'*X';
%
Q0 = X*A+A'*X';
L = eye(n,n);
N = B2;
R = C2;
%
LMImanual = [Q0+L*X*N*Y0*R+L*X0*N*Y*R-L*X0*N*Y0*R+(L*X*N*Y0*R+L*X0*N*Y*R-L*X0*N*Y0*R)',...
             L*(X-X0)*N+R'*(Y-Y0)'*G';...
            (L*(X-X0)*N+R'*(Y-Y0)'*G')',...
             -(G+G')];
         
% % BMI 制約式の左辺
disp("# 制約がブロック行列で表されない")

% Fstr = "X*A+X*B2*Y*C2+A'*X'+C2'*Y'*B2'*X'";
Fstr = "-X*(-A-B2*Y*C2)+(A+(-B2)*(-Y)*C2)'*X'";
% Fstr = "eye(n)*X*(A+B2*Y*C2)+(A+B2*Y*C2)'*X'*eye(n,n)+(zeros(n)-zeros(n))*(zeros(n)+zeros(n))";

tStart = tic;
[LMIauto, BMIauto] = bmiparser(Fstr,{'X','Y'},{'X0','Y0'});
tEnd = toc(tStart)

%
BMItt = BMImanual - BMIauto;
LMItt = LMImanual - LMIauto;

disp("BMItt: "+evaltest(BMItt));
disp("LMItt: "+evaltest(LMItt));

cll = cat(1,cll,{"9. 極配置",evaltest(BMItt),evaltest(LMItt),tEnd});


%% パターンテストの結果
disp(newline)
TestTable = cell2table(cll,'VariableNames',{'Pattern','BMItt','LMItt','RunTime'})



function eval = evaltest(tt)
    if isequal(tt,zeros(size(tt)))
        eval = true;
    else
        eval = false;
    end
end




