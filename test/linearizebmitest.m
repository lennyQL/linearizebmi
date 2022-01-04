%% 出力フィードバック制御のBMI制約についてのテスト
%   linearizebmiを評価するテストプログラム
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

% Hinf問題の制約の左辺，linearizebmiのための制約の記述
Fstr = "[X*A+X*B2*Y*C2+A'*X'+C2'*Y'*B2'*X'  X*B1+X*B2*Y*D21  C1'+C2'*Y'*D12';"+...
        "B1'*X'+D21'*Y'*B2'*X'               -I               D11'+D21'*Y'*D12';"+...
        "C1+D12*Y*C2                         D11+D12*Y*D21    -I]";

% parser(自作した関数)
tStart = tic;
[LMIauto, ~,~, BMIauto] = linearizebmi(Fstr,{'X','Y'},{'X0','Y0'});
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
[LMIauto, ~,~, BMIauto] = linearizebmi(Fstr,{'X','Y'},{'X0','Y0'});
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
[LMIauto, LMIstr, g, BMIauto] = linearizebmi(Fstr,{'X','Y'},{'X0','Y0'});
tEnd = toc(tStart)

LMIauto
LMIstr
% g
% g.sdpvar
% g.data
% g.str
% g.str.Q
% g.str.L

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
[LMIauto, ~,~, BMIauto] = linearizebmi(Fstr,{'X','Y'},{'X0','Y0'});
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
[LMIauto, ~,~, BMIauto] = linearizebmi(Fstr,{'X','Y'},{'X0','Y0'});
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
[LMIauto, ~,~, BMIauto] = linearizebmi(Fstr,{'X','Y'},{'X0','Y0'});
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
[LMIauto, ~,~, BMIauto] = linearizebmi(Fstr,{'X','Y'},{'X0','Y0'});
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
[LMIauto, ~,~, BMIauto] = linearizebmi(Fstr,{'X','Y'},{'X0','Y0'});
tEnd = toc(tStart)

BMItt = BMImanual - BMIauto;
LMItt = LMImanual - LMIauto;

disp("BMItt: "+evaltest(BMItt));
disp("LMItt: "+evaltest(LMItt));

cll = cat(1,cll,{"8. blkbmiによる表現",evaltest(BMItt),evaltest(LMItt),tEnd});



%% パターン9 分割行列Gも決定変数のとき（決定変数3つ）
disp(newline)
disp("###*** パターン9 ***###")
disp("# 分割行列も決定変数のとき")

Z = sdpvar(size(Y,1),size(Y,1),'full');
Z0 = rand(size(Z));

Fstr = "[X*(A+B2*Y*C2)+(A+B2*Y*C2)'*X'  X*(B1+B2*Y*D21)     (C1+D12*Y*C2)';"+...
        "(B1+B2*Y*D21)'*X'               -eye(p1)           (D11+D12*Y*D21)';"+...
        "C1+D12*Y*C2                     D11+D12*Y*D21       -eye(p1)]";

% [LMIauto, Lstr,~, BMIauto] = linearizebmi(Fstr,{'X','Y'},{'X0','Y0'},'G')
tStart = tic;
[LMIauto, Lstr,~, BMIauto] = linearizebmi(Fstr,{'X','Y','Z'},{'X0','Y0','Z0'},'G')
tEnd = toc(tStart)


LMImanualZ = [Q0+L*X*N*Y0*R+L*X0*N*Y*R-L*X0*N*Y0*R+...
    (L*X*N*Y0*R+L*X0*N*Y*R-L*X0*N*Y0*R)',...% (1,1)
     L*(X-X0)*N+R'*(Y-Y0)'*Z0',...          % (1,2)
    (G*(Y-Y0)*R)';...                       % (1,3)
    (L*(X-X0)*N+R'*(Y-Y0)'*Z0')',...        % (2,1)
     -(Z+Z'),...                            % (2,2)
     Z-Z0;...                               % (2,3)
     G*(Y-Y0)*R,...                         % (3,1)
    (Z-Z0)',...                             % (3,2)
     -(G+G')];                              % (3,3)

BMItt = BMImanual - BMIauto;
LMItt = LMImanualZ - LMIauto;

disp("BMItt: "+evaltest(BMItt));
disp("LMItt: "+evaltest(LMItt));



%% 極配置問題, 制約がブロック行列で表されない場合
disp(newline)
disp("<<<======***** 極配置問題 *****======>>>")


X=sdpvar(n,n);
Y=sdpvar(m2,p2, 'full');

G = eye(size(Y,1));


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
[LMIauto, LMIstr,~, BMIauto] = linearizebmi(Fstr,{'X','Y'},{'X0','Y0'});
% G=sdpvar(size(Y,1),size(Y,1),'full');
% H=rand(size(G));
% [LMIauto, Lstr,~, BMIauto] = linearizebmi(Fstr,{'X','Y','G'},{'X0','Y0','H'})
tEnd = toc(tStart)


%
LMIstr
BMItt = BMImanual - BMIauto
LMItt = LMImanual - LMIauto

disp("BMItt: "+evaltest(BMItt));
disp("LMItt: "+evaltest(LMItt));

cll = cat(1,cll,{"9. 極配置",evaltest(BMItt),evaltest(LMItt),tEnd});

%% 構文errorテスト
disp(newline)
disp("###*** パターン ***###")
disp("# 構文errorテスト")

X=sdpvar(n,n);
Y=sdpvar(m2,p2, 'full');

try 
    X0=sdpvar(n,n);
    Y0=sdpvar(m2,p2, 'full');
    Fstr = "[X*(A+B2*Y*C2)+(A+B2*Y*C2)'*X';";
    [LMIauto, ~,~, BMIauto] = linearizebmi(Fstr,{'X','Y'},{'X0','Y0'});
    disp('correct')
catch ME
    disp([ME.identifier ME.message]);
end

Fstr = "X*(A+B2*Y*C2)+(A+B2*Y*C2)'*X'"; 
try 
    X0=sdpvar(n,p2);
    Y0=sdpvar(m2,p2, 'full');
    LMIauto = linearizebmi(Fstr,{'X','Y'},{'X0','Y0'});
    disp('correct')
catch ME
    disp(ME.message);
end
try 
    X0=sdpvar(n,n);
    Y0=sdpvar(n,n, 'full');
    LMIauto = linearizebmi(Fstr,{'X','Y'},{'X0','Y0'});
    disp('correct')
catch ME
    disp(ME.message);
end

try 
    X0=sdpvar(n,n);
    Y0=sdpvar(m2,p2, 'full');
    LMIauto = linearizebmi(Fstr,{'X','Y'},{'X0','Y0'},G);
    disp('correct')
catch ME
    disp([ME.identifier ME.message]);
end

try 
    X=eye(n);
    Y=sdpvar(m2,p2, 'full');
    LMIauto = linearizebmi(Fstr,{'X','Y'},{'X0','Y0'});
    disp('correct')
catch ME
    disp([ME.identifier ME.message]);
end

try 
    X=sdpvar(n,n);
    Y=sdpvar(m2,p2, 'full');
    Z=sdpvar(m2,p2);
    linearizebmi(Fstr,{'X','Z'},{'X0','Y0'});
    disp('correct')
catch ME
    disp([ME.identifier ME.message]);
end

%% デフォルト
disp(newline)
disp("###*** パターン ***###")
disp("# デフォルト")

X0=sdpvar(n,n);
Y0=sdpvar(m2,p2, 'full');

Fstr = "[X*(A+B2*Y*C2)+(A+B2*Y*C2)'*X'  X*(B1+B2*Y*D21)     (C1+D12*Y*C2)';"+...
        "(B1+B2*Y*D21)'*X'               -eye(p1)           (D11+D12*Y*D21)';"+...
        "C1+D12*Y*C2                     D11+D12*Y*D21       -eye(p1)]";

% LMIauto = linearizebmi(Fstr,{'X','Y'},{'X0','Y0'})
[LMIauto, LMIstr, out] = linearizebmi(Fstr,{'X','Y'},{'X0','Y0'},'G')

%% そもそもLMIの場合
disp(newline)
disp("###*** パターン ***###")
disp("# LMI")

X=sdpvar(n,n);
Y=sdpvar(m2,p2, 'full');
X0=rand(size(X));
Y0=rand(size(Y));

g = sdpvar;
h = sdpvar;
g0 = rand(size(g));
h0 = rand(size(h));

% diag(g)
% dig = diag(X)


% Fstr = "A*X+(A*X)'+A*X+B2*Y*B2'+(B2*Y*B2')'";
% Fstr = "-X";
Fstr = "trace(X)*trace(Y)+(trace(X)*trace(Y))'"
% Fstr = "X*(A+B2*Y*C2)+(A+B2*Y*C2)'*X'"

% Fstr = "diag(g)*h+(diag(g)*h)'"

[LMIauto, ~, out] = linearizebmi(Fstr,{'trace(X)','trace(Y)'},{'trace(X0)','trace(Y0)'})
% [LMIauto, ~, out] = linearizebmi(Fstr,{'diag(g)','h'},{'diag(g0)','h0'})
% [LMIauto, ~, out] = linearizebmi(Fstr,{'X','Y'},{'X0','Y0'})

is(LMIauto,'linear')

out.sdpvarname


%% Webマニュアル用 (極配置)
disp(newline)
disp("###*** Webマニュアル用 (極配置) ***###")

P=sdpvar(n,n);
K=sdpvar(m2,p2, 'full');
P0=rand(size(X));
K0=rand(size(Y));
% P0=eye(size(X));
% K0=eye(size(Y));



% Fstr = "X*A+X*B2*Y*C2+A'*X'+C2'*Y'*B2'*X'";
Fstr = "P*(A+B2*K*C2)+(P*(A+B2*K*C2))'";

G=eye(size(G));
[LMIauto, LMIstr] = linearizebmi(Fstr,{'P','K'},{'P0','K0'},'G')

G=sdpvar(size(Y,1),size(Y,1),'full');
H=rand(size(G));
[LMIauto, Lstr] = linearizebmi(Fstr,{'P','K','G'},{'P0','K0','H'})

G=sdpvar(size(Y,1),size(Y,1),'full');
H=rand(size(G));
M=sdpvar(m2,m2);
M0=eye(size(M));
[LMIauto, Lstr] = linearizebmi(Fstr,{'P','K','G','M'},{'P0','K0','H','M0'})

% LMIstr
% Lstr


%% Input数拡張(オプション付き)
disp(newline)
disp("###*** method選択(オプション付き) ***###")

P=sdpvar(n,n);
K=sdpvar(m2,p2, 'full');
P0=rand(size(P));
K0=rand(size(K));

opts = linearizebmiOptions;


% Fstr = "P*(A+B2*K*C2)+(P*(A+B2*K*C2))'";
Fstr = "[P*(A+B2*K*C2)+(A+B2*K*C2)'*P'  P*(B1+B2*K*D21)     (C1+D12*K*C2)';"+...
        "(B1+B2*K*D21)'*P'               -eye(p1)           (D11+D12*K*D21)';"+...
        "C1+D12*K*C2                     D11+D12*K*D21       -eye(p1)]";


% method:0
disp("-- method:0 --")
opts.method = 0;
[LMIauto, LMIstr] = linearizebmi(Fstr,{'P','K'},{'P0','K0'},'',opts)

% method:1
disp("-- method:1 --")
G=sdpvar(m2,m2,'full');
G0=eye(size(G));
opts.method = 1;
[LMIauto, LMIstr] = linearizebmi(Fstr,{'P','K','G'},{'P0','K0','G0'},'',opts)

% method:2
disp("-- method:2 --")
G=eye(size(G));
opts.method = 2;
[LMIauto, LMIstr] = linearizebmi(Fstr,{'P','K'},{'P0','K0'},'G',opts)

% method:3
disp("-- method:3 --")
G=sdpvar(m2,m2);
G0=eye(size(G));
opts = linearizebmiOptions('method',3);
[LMIauto, LMIstr] = linearizebmi(Fstr,{'P','K','G'},{'P0','K0','G0'},'',opts)

% method:4
disp("-- method:4 --")
G=sdpvar(m2,m2,'full');
G0=eye(size(G));
M=sdpvar(m2,m2);
M0=eye(size(M));
opts = linearizebmiOptions(opts,'method',4);
[LMIauto, LMIstr] = linearizebmi(Fstr,{'P','K','G','M'},{'P0','K0','G0','M0'},'',opts)



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




