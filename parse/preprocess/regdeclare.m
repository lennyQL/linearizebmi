%% 正規表現用マジックナンバーの宣言(初期化)する関数，正規表現を要素ごとに分解
function regdeclare()
    clear regdeclare

    % 変数名，転置付き
    % ex) A, B1'
    global VAR
    VAR = "\w+(\')*";
    
    % 関数
    % ex) eye(n,n)
    global FUNC
    % FUNC = "\w+\([\w\,\'\-\+\*]+\)(\')*";
    FUNC = "\w+\((?:[^\(\)]|\w+\([^\)]*\))*\)(\')*";
    
    % 関数もしくは変数名, 値を持つトークン
    % ex) eye(n), A'
    global FUNC_OR_VAR
    FUNC_OR_VAR = "(\-)*"+"("+FUNC+"|"+VAR+")";
    
    % 多項式の項
    % ex) A, A*B1, eye(n)*A'
    global TERM
    TERM = FUNC_OR_VAR+"(\*"+FUNC_OR_VAR+")*";
    
    % 多項式
    % ex) A+B1'*eye(n)
    global POLY
    POLY = FUNC_OR_VAR+"([\+\-\*]"+FUNC_OR_VAR+")*";
    
    
    % 括弧の前後に積がない & 関数の括弧でない
    % ex) (A+B*C)
    % eye(n)とはマッチングしない
    global BRACKET_POLY
    BRACKET_POLY = "(?<!\w+)"+"\("+POLY+"\)";
    
    % 括弧の転置
    % ex) (A+B*C)'
    global BRACKET_TRANSPOSE
    BRACKET_TRANSPOSE = BRACKET_POLY+"\'";
    
    % 括弧と括弧の積
    % ex) (A+B*C)*(A+B*C)
    % eye(n)*(A+B*C) とはマッチングしない
    global PROD_BRACKET 
    global BRACKET_PROD_BRACKET
    PROD_BRACKET = "(\*"+BRACKET_POLY+")+"+"(?!.*\))";
    BRACKET_PROD_BRACKET = BRACKET_POLY + "(\*"+BRACKET_POLY+")+";
    
    % 括弧の前に積
    % ex) D*(A+B*C) 
    global LEFT_PROD 
    global LEFT_PROD_BRACKET
    LEFT_PROD = "("+FUNC_OR_VAR+"\*)+";
    LEFT_PROD_BRACKET = LEFT_PROD + BRACKET_POLY;
    
    % 括弧の後に積
    % ex) (A+B*C)*D
    % eye(n)*D とはマッチングしない
    global PROD_RIGHT 
    global BRACKET_PROD_RIGHT
    PROD_RIGHT = "(\*"+FUNC_OR_VAR+")+"+"(?!.*\))";
    BRACKET_PROD_RIGHT = BRACKET_POLY + "(\*"+FUNC_OR_VAR+")+";
    
    
    % 縦ベクトルと横ベクトルの積
    global VEC
    global PROD_VEC
    global VEC_PROD_VEC
    VEC = "\[[^\[]*?\]";
    PROD_VEC = "\*" + VEC;
    VEC_PROD_VEC = VEC + PROD_VEC; 
    
    
    % 括弧の多項式
    % ex) (A+B*C)*A, A*B, A+(B+C*K)
    global BRACKET_POLY_OR_POLY
    global POLY_BRACKET
    BRACKET_POLY_OR_POLY = "("+BRACKET_POLY+"|"+POLY+")";
    POLY_BRACKET = BRACKET_POLY_OR_POLY+"((+|-|*)"+BRACKET_POLY_OR_POLY+")*";
    
    % 行列の転置
    global MAT_TRANSPOSE
    global MAT_COL
    MAT_TRANSPOSE = "("+ "(?<!\w+)"+"\("+VEC+"\)"+"\'" +"|" +VEC+"\'" +")";
    % MAT_COL = "("+BRACKET_POLY_OR_POLY+"(\s)*)+"+"(?=(;|]))";
    MAT_COL = "(?<=([|;))"+"[^;]+"+"(?=(;|]))";
    
end