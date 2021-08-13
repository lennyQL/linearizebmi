function [LMI,BMI,Q0,L,N,R] = bmiparser(S, vlist, v0list, G)
% BMIの文字列を受け取り，逐次LMIの値を返すパーサー
%
%   S: BMIの文字列
%       ex) "P*A+P*B*K*C+A'*P'+C'*K'*B'*P'"
%   vlist: 決定変数の文字列
%       ex) ["P","K"]
%   v0list: 暫定解の文字列
%       ex) ["P0","K0"]
%   G: gammaの定数倍，デフォルトで単位行列
%
%
%   ※ワークスペースの変数の値を利用して計算するので，
%     bmiparser関数を呼び出す前に，
%     あらかじめに変数を宣言する必要がある．
%   ex) 
%       X = sdpvar(2,2)
%       Y = sdpvar(3,2)
%       A = rand(2,3)
%       X0= rand(2,2) 
%       Y0= rand(3,2)
%       LMI = bmiparser("X*A*Y+Y'*A'*X'",["X","Y"],["X0","Y0"])
%
%
%   現段階：
%       ・演算子は，+,*,' ,しか実装できていない
%       ・()は使えない
%       ・スカラー倍にまだ対応できていない
%

% 文字列を文字ベクトルに変換
if isa(S,'string')
    S = char(S);
end
% str2sym(S)

% 関数引数の取得
Xstr =vlist{1};
Ystr =vlist{2};
X0str=v0list{1};
Y0str=v0list{2};


%% 字句解析
% ex) "P*A+P*B*K*C"
% =>  {{"P","B"},{"P","B","K","C"}}
% ex) "P'*A'"
% =>  {{"P'","A'"}}

termlist = {}; % 項のリスト
varlist = {};  % 変数定数の配列(一時的)
varstr = "";   % 1つの変数の文字列(一時的)

% 制約が行列の場合
columlist = {}; % 行のリスト，termlistをこの中に入れる
colnum = 1;     % 行数
rownum = 1;     % 列数

% デバッグ用:
% disp("==========")
% 入力文字列を1文字ずつ解析
for i=1:strlength(S)
    % デバッグ用:
    % disp("-----")
    % S(i),varstr,varlist,termlist

    
    % 変数の文字列の場合
    % ex) 'var' = 'v'+'a'+'r' 
    % 正規表現[a-zA-Z_0-9]
    % ex) a0は変数だが，0aは変数でない -> どう処理する？
    if regexp(S(i), '\w')
        % どこまでが変数名かを推定する
        % 変数名の更新
        varstr = varstr + S(i);
        if i == strlength(S) 
            % 終端文字の場合，varlistをtermlistに追加
            % "A]"とかだと別処理が必要，"[...]"の処理
            [varlist, varstr] = updateList(varlist,varstr);
            [termlist, varlist] = updateList(termlist,varlist,1);
            [columlist, termlist] = updateList(columlist,termlist);
        end
        continue
    end
    
    % 演算子の場合
    % 和+
    if regexp(S(i), '\+') 
        % varstrをvarlistに追加, varstr初期化
        [varlist, varstr] = updateList(varlist,varstr);
        % varlistをtermlistに追加，varlist初期化
        [termlist, varlist] = updateList(termlist,varlist,1);
        continue
    end
    % 積*
    if regexp(S(i), '\*')    
        [varlist, varstr] = updateList(varlist,varstr);
        continue
    end
    
    % 差-
    % ex) "-A*B": {"-","A","B"}
    if regexp(S(i), '\-')
        % 前までの項をリスト追加(あれば)
        if varstr ~= ""
            [varlist, varstr] = updateList(varlist,varstr);
            [termlist, varlist] = updateList(termlist,varlist,1);
        end
        % 新しい項の1文字目に入れる
        varstr = varstr + S(i);
        [varlist, varstr] = updateList(varlist,varstr);
        continue
    end
    
    % シングルクオーテーションの場合
    % ex) A'*B'：行列の転置をする時に使う
    % 現状 \w と同じ処理
    if regexp(S(i), "\'")
        varstr = varstr + S(i);
        if i == strlength(S) 
            % 終端文字の場合，termstrをtermlistに追加
            [varlist, varstr] = updateList(varlist,varstr);
            [termlist, varlist] = updateList(termlist,varlist,1);
            [columlist, termlist] = updateList(columlist,termlist);
        end
        continue
    end
        
    % 空文字の場合
    % ex) '[A B]' : 行列の列を生成するときに使う
    if regexp(S(i), '\s')    
        if ~isempty(regexp(S(i-1), '\s', 'once')) ||...
           ~isempty(regexp(S(i), '[\*\+]', 'once'))   
            continue
        elseif regexp(S(i-1), "[\w\']")
            [varlist, varstr] = updateList(varlist,varstr);
            [termlist, varlist] = updateList(termlist,varlist,1);
            [columlist, termlist] = updateList(columlist,termlist);
            % 列数の更新
            rownum = rownum + 1;
            continue
        end
    end
    % セミコロン;の場合
    % ex) '[A;B]' : 行列の行を生成するときに使う
    if regexp(S(i), ';')
        [varlist, varstr] = updateList(varlist,varstr);
        [termlist, varlist] = updateList(termlist,varlist,1);
        [columlist, termlist] = updateList(columlist,termlist);
        % 行数更新
        colnum = colnum + 1;
        rownum = 1;
        continue
    end
    
    % [の場合
    % 制約は行列
    if regexp(S(i), '\[')
        % disp("Helo")
    end
    % ]の場合
    % 行列おわり
    if regexp(S(i), '\]')
        % disp("bye")
        [varlist, varstr] = updateList(varlist,varstr);
        [termlist, varlist] = updateList(termlist,varlist,1);
        [columlist, termlist] = updateList(columlist,termlist);
    end
end
% smatrixの生成
% Sの最終形態(cellの行列), columlistを縦に並べたもの
if colnum > 1 || rownum > 1
    % ブロック行列の場合
    smatrix = reshape(columlist,colnum,rownum).';
else
    % そうでない場合
    smatrix = columlist;
end

% columlist
% smatrix = cat(2,cat(2,smatrix{1,1},smatrix{2,1}),smatrix{3,1});
% smatrix = reshape(smatrix,[colnum,rownum])

% デバッグ用:
% smatrix
% smatrix{1,1}
% smatrix{1,1}{1,1}
% smatrix{1,1}{1,1}{1,1}
% smatrix{2,1}{1,2}
% smatrix{2,1}{1,2}{1,1}
% smatrix{2,1}{1,2}{1,2}
% smatrix{3,1}{1,3}{1,2}

% smatrix(1,2)
% smatrix{1,2}
% columlist{1,3}{2,1}
% smatrix{1,3}{2,1}


%% heの分離

orgmatrix = smatrix; % 転置なし行列
hematrix = smatrix;  % 転置あり行列

orgtermlist = {};    % 項のリスト(初期化)
hetermlist = {};     % 項の転置ありリスト(初期化)

for col=1:size(smatrix,1)
    % 各行ベクトル
    for row=1:size(smatrix,2)
        % 各行列要素
        % disp(col+" "+row)
        termlist = smatrix{col,row};
        for i=1:size(termlist,1)
            % 要素の各項
            term = termlist{i,1};
            for j=1:size(term,2)
                % 項のそれぞれの変数
                var = term{1,j};
                if regexp(var, "\'")
                    % var
                    % hetermlist: 転置あり
                    hetermlist = updateList(hetermlist,term,1);
                    break
                else
                    % termlist: 転置なし
                    % if col == row && col >= 2
                    %     % (2,2),(3,3)要素は1/2倍する
                    %     % スカラー倍を実装しないとできない
                    % end
                    orgtermlist = updateList(orgtermlist,term,1);
                    break
                end
            end
        end
        % 転置なし行列の更新
        orgmatrix(col,row) = {orgtermlist};
        orgtermlist = {};
        % 転置あり行列の更新
        hematrix(col,row) = {hetermlist};
        hetermlist = {};
    end
end
% 'があったら，hematrixに追加
% デバッグ用:
% orgmatrix
% hematrix


%% 双線形項の分離
% 行列
linearmatrix = orgmatrix; % 定数項，1次項
binearmatrix = orgmatrix; % 双線形項

% 一時的リスト
linearlist = {}; % 定数項，1次項
binearlist = {}; % 双線形項


for col=1:size(orgmatrix,1)
    % 各行ベクトル
    for row=1:size(orgmatrix,2)
        % 各行列要素
        % disp(col+" "+row)
        termlist = orgmatrix{col,row};
        for i=1:size(termlist,1)
            % 要素の各項
            term = termlist{i,1};
            %
            varcount = 0; % 決定変数の数
            for j=1:size(term,2)
                % 項のそれぞれの変数
                var = term{1,j};
                if ~isempty(regexp(var, Xstr, 'once')) ||...
                   ~isempty(regexp(var, Ystr, 'once'))
                    varcount = varcount + 1;
                end
                if varcount >= 2
                    binearlist = updateList(binearlist,term,1);
                    break
                elseif j == size(term,2)
                    linearlist = updateList(linearlist,term,1);
                    break
                end
            end
        end
        % 線形項の行列
        linearmatrix(col,row) = {linearlist};
        linearlist = {};
        % 双線形項の行列
        binearmatrix(col,row) = {binearlist};
        binearlist = {};
    end
end
% デバッグ用:
% linearmatrix
% binearmatrix

% binearmatrix{1,1}{1,1}
% binearmatrix{1,2}{1,1}


%% BMI一般化，Q,L,N,Rの取得
Q = linearmatrix; % 定数項，一次項
L = {}; % 双線形項の定数行列(左)
N = {}; % 双線形項の定数行列(中)
R = {}; % 双線形項の定数行列(右)


for col=1:size(binearmatrix,1)
    % 各行ベクトル
    for row=1:size(binearmatrix,2)
        % 各行列要素
        % disp(col+" "+row)
        termlist = binearmatrix{col,row};
        %
        xidx = 0; % 双線形項におけるXstrの位置
        yidx = 0; % 双線形項におけるYstrの位置
%         Llist = {}; % 一時リスト
%         Nlist = {}; 
%         Rlist = {}; 
        for i=1:size(termlist,1)
            % 要素の各項
            term = termlist{i,1};
            for j=1:size(term,2)
                % 項のそれぞれの変数
                var = term{1,j};
                if regexp(var, Xstr, 'once')
                    xidx = j;
                elseif regexp(var, Ystr, 'once') 
                    yidx = j;
                end
            end    
            % P,K(X,Y)でリストを3分割する
            l = term(1:xidx-1);
            n = term(xidx+1:yidx-1);
            r = term(yidx+1:end);
            % "PK"の場合： "1*P*1*K*1"と処理する
            if isempty(l) 
                l = {"1eye"};
            elseif l{1,1} == "-"
                l = ["-", "1eye"];
            end
            if isempty(n)
                n = {"1eye"};
            elseif n{1,1} == "-"
                n = ["-", "1eye"];
            end
            if isempty(r)
                r = {"1eye"};
            elseif r{1,1} == "-"
                r = ["-", "1eye"];
            end

            % L,N,Rの更新
            % 未完成(仮)
            % 双線形項が1行目にある場合にしか対応できない
            if col == 1                
                if isempty(L)
                    L = updateList(L,l);
                end
                N = {n};
                R = updateList(R,r,1);
            end
        end
        
        % 未完成(仮)
        if col == 1 && isempty(termlist)
            R = updateList(R,{"0zero"},1);
        end
        
%         L = updateList(L,Llist);
%         N = updateList(N,Nlist);
%         R = updateList(R,Rlist);
    end
    
    % 未完成(仮)
    if col >= 2
        L = updateList(L,{"0zero"},1);
    end
end
% デバッグ用:
% L{1,1}
% isempty(L{1,1})
% Q,Q{3,1}{2,1}
% L,N,R


%% 逐次LMIに変形してCalc
% evalinを使う，workspaceの変数の値の取得


% 決定変数の取得
X = evalin('base', Xstr);
Y = evalin('base', Ystr);
% 暫定解の取得
X0 = evalin('base', X0str);
Y0 = evalin('base', Y0str);


% ブロック行列のサイズの計算
% 制約行列の対角ブロックから割り出す
colsize = []; % 各ブロックの行サイズ
rowsize = []; % 各ブロックの列サイズ
for i=1:size(orgmatrix,1)
    idx = 1;
    s = orgmatrix{i,i}{1,1}{1,1};
    for j=1:size(orgmatrix,1)
        if s == "-"
            idx = idx+1;
            break
        end
    end
    s = orgmatrix{i,i}{1,1}{1,idx}; % 項の一番左の変数
    e = orgmatrix{i,i}{1,1}{1,end}; % 項の一番右の変数
    ssize = size(evalin('base',s),1);
    esize = size(evalin('base',e),2);
    colsize = cat(2,colsize,ssize);
    rowsize = cat(2,rowsize,esize);
end
% colsize,rowsize


% 線形項の計算
Qeval = [];
for col=1:size(Q,1)
    % 各行ベクトル
    Qcol = [];
    for row=1:size(Q,2)
        % 各行列要素
        % disp(col+" "+row)
        termlist = Q{col,row};
        %
        Qevalelem = 0;
        for i=1:size(termlist,1)
            term = termlist{i,1};
            qeval = 1;
            for j=1:size(term,2)
                var = term{1,j};
                if var == "-"
                    qeval = -qeval;
                else
                    qeval = qeval * evalin('base', var);
                end
            end
            Qevalelem = Qevalelem + qeval;
        end
        
        % 要素に項がない場合，ゼロ行列を入れる
        if isempty(termlist)
            % col,row
            z = zeros(colsize(col),rowsize(row));
            Qcol = cat(2,Qcol,z);
        else
            % 対角要素は1/2倍
            if col == row && col >= 2
                % (1,1)要素に対しても実行する必要がある
                % 未実装
                Qevalelem = Qevalelem / 2;
            end
            Qcol = cat(2,Qcol,Qevalelem);
        end
        % Qcol
        
    end
    Qeval = cat(1,Qeval,Qcol);
end
% Qeval


% 双線形項の左定数Lの計算
Leval = [];
for i=1:size(L,1)
    term = L{i,1};
    leval = 1; 
    for j=1:size(term,2)
        var = term{1,j};
        if var == "1eye"
            leval = leval * eye(size(X,1));
        elseif var == "0zero"
            leval = leval * zeros(rowsize(i),size(Leval,2));
        elseif var == "-"
            leval = -leval;
        else 
            leval = leval * evalin('base', var);
        end
    end
    % 列(縦)ベクトルを生成
    Leval = cat(1,Leval,leval);
end
% Leval


% 双線形項の中定数Nの計算
Neval = [];
for i=1:size(N,1)
    term = N{i,1};
    neval = 1; 
    for j=1:size(term,2)
        var = term{1,j};
        if var == "1eye"
            neval = neval * eye(size(X,2),size(Y,1));
        elseif var == "-"
            neval = -neval;
        else
            neval = neval * evalin('base', var);
        end
    end
    Neval = neval;
end
% Neval


% 双線形項の右定数Rの計算
Reval = [];
for i=1:size(R,1)
    term = R{i,1};
    reval = 1; 
    for j=1:size(term,2)
        var = term{1,j};
        if var == "1eye"
            reval = reval * eye(size(Y,2));
        elseif var == "0zero"
            reval = reval * zeros(size(Reval,1),rowsize(i));
        elseif var == "-"
            reval = -reval;
        else
            reval = reval * evalin('base', var);
        end
    end
    % 行(横)ベクトルを生成
    Reval = cat(2,Reval,reval);
end
% Reval


% 双線形項の計算, LXNYR
Bieval = Leval * X * Neval * Y * Reval;

% BMIの値の計算, 多分使わない，デバッグ用
BMIeval = Qeval + Bieval;


if nargin<4
    G=eye(size(Y,1));
elseif isa(G,'string')
    G=evalin('base',G);
end

% 拡大した LMI 条件, heなし
LMIeval=[Qeval+replace(Bieval,Y,Y0)+replace(Bieval,X,X0)-Leval*X0*Neval*Y0*Reval,...
         Leval*(X-X0)*Neval;...
         G*(Y-Y0)*Reval,...
         -G];
% heあり
% LMIeval = [Qeval+Leval*X*Neval*Y0*Reval+Leval*X0*Neval*Y*Reval-Leval*X0*Neval*Y0*Reval+...
%         (Qeval+Leval*X*Neval*Y0*Reval+Leval*X0*Neval*Y*Reval-Leval*X0*Neval*Y0*Reval)',...
%          Leval*(X-X0)*Neval+Reval'*(Y-Y0)'*G';...
%         (Leval*(X-X0)*Neval+Reval'*(Y-Y0)'*G')',...
%          -2*G];


     
%% デバッグ用出力
%%%% heなし
Q0= Qeval; % 線形項
L = Leval; % 双線形項の定数行列(左)
N = Neval; % 双線形項の定数行列(中)
R = Reval; % 双線形項の定数行列(右)


%% 関数の出力
LMI = LMIeval + LMIeval';
BMI = BMIeval + BMIeval';



end



%% 内部関数群
% Listの更新（追加）と初期化
function [L,V] = updateList(L,V,n)
    % ex) 
    %   LにVを追加，Vを初期化:
    %     [L,V] = updateList(L,V)
    %   LにVを追加のみ: 
    %     L = updateList(L,V)

    if nargin<3
        n = 2;
    end

    % Listの追加
    % n=1: 縦に追加
    % n=2: 横に追加
    L = cat(n,L,{V});
    
    if isa(V, "string")
        % stringの初期化
        V = "";
    elseif isa(V, "cell")
        % cellの初期化
        V = {};
    end
end

% Matrixの更新
% function [M,L] = updateMatrix(M,L,col,row)
%     collist = M{col,1};
%     rowlist = M{col,1}{1,row};
%     
%     rowlist = updateList(rowlist,L);
%     collist{1,row} = rowlist;
%     M{col,1} = collist;
%     
%     L = {};
% end



