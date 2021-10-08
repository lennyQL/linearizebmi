%% 括弧()を展開する関数
function S = divbracket(S)

    % 対応する正規表現
    % TODO:
    %   ブロック行列を括弧で囲んだ表現 ([...])'をどう処理するか
    %   変数の前にマイナス-があるときどう処理するか

    
    % 正規表現用変数(global)を呼び出す
%     regdeclare();
    global FUNC_OR_VAR
    global TERM
    global POLY
    global BRACKET_POLY
    global BRACKET_TRANSPOSE
    global PROD_BRACKET 
    global BRACKET_PROD_BRACKET
    global LEFT_PROD 
    global LEFT_PROD_BRACKET
    global PROD_RIGHT 
    global BRACKET_PROD_RIGHT
    
    
    % 括弧の転置，括弧内の積の順序が逆転
    if regexp(S,BRACKET_TRANSPOSE)
        % 正規表現に対応する文字列
        bracket = regexp(S,BRACKET_TRANSPOSE,'match','once');
        % 括弧の中の式
        bracketin = regexp(bracket,POLY,'match','once');
        %
        % 括弧の中の各項
        term = regexp(bracketin,TERM,'match');
        % 括弧の中の文字列(最終的に求めたいもの)
        bracketstr = "";
        for i=1:length(term)
            % 各項の文字列
            termstr = "";
            % 各項の変数のcell配列
            var = regexp(term{1,i},FUNC_OR_VAR,'match');
            for j=1:length(var)
                % 配列を逆から処理，転置'を加えて*で連結
                if j==1
                    termstr = termstr + string(var(end-j+1)) + "'";
                else
                    termstr = termstr + "*" + string(var(end-j+1)) + "'";
                end
            end

            % 上で処理した各項を+で連結
            if i==1
                bracketstr = bracketstr + termstr;
            else
                bracketstr = bracketstr + "+" + termstr;
            end
        end
        % 括弧()で囲み置換する
        rep = "(" + bracketstr + ")";
        S = regexprep(S,BRACKET_TRANSPOSE,rep,'once');
    
    
    % 括弧と括弧の積
    elseif regexp(S,BRACKET_PROD_BRACKET)
        % 正規表現と対応する文字列
        bracketterm = regexp(S,BRACKET_PROD_BRACKET,'match','once');
        % 括弧との積
        mult = regexp(bracketterm,PROD_BRACKET,'match','once');
        % 括弧の項
        bracket = regexp(bracketterm,BRACKET_POLY,'match','once');
        % 括弧の中の式
        bracketin = regexp(bracket,POLY,'match','once');
        % 置き換える文字列，括弧の中の各項との積
        rep = regexprep(bracketin,TERM,"$0"+mult);
        % 括弧を展開した式に変換する
        S = regexprep(S,BRACKET_PROD_BRACKET,rep,'once');
        
        
    % 括弧の前に積
    elseif regexp(S,LEFT_PROD_BRACKET)
        % 正規表現と対応する文字列
        bracketterm = regexp(S,LEFT_PROD_BRACKET,'match','once');
        % 括弧との積
        mult = regexp(bracketterm,LEFT_PROD,'match','once');
        % 括弧の項
        bracket = regexp(bracketterm,BRACKET_POLY,'match','once');
        % 括弧の中の式
        bracketin = regexp(bracket,POLY,'match','once');
        % 置き換える文字列，括弧の中の各項との積
        rep = regexprep(bracketin,TERM,mult+"$0");
        % 括弧を展開した式に変換する
        S = regexprep(S,LEFT_PROD_BRACKET,"("+rep+")",'once');
        
        
    % 括弧の後に積
    elseif regexp(S,BRACKET_PROD_RIGHT)
        bracketterm = regexp(S,BRACKET_PROD_RIGHT,'match','once');
        mult = regexp(bracketterm,PROD_RIGHT,'match','once');
        bracket = regexp(bracketterm,BRACKET_POLY,'match','once');
        bracketin = regexp(bracket,POLY,'match','once');
        rep = regexprep(bracketin,TERM,"$0"+mult);
        S = regexprep(S,BRACKET_PROD_RIGHT,rep,'once');
    
        
    % 括弧の前後に積がない，括弧を外す
    elseif regexp(S,BRACKET_POLY)
        bracketterm = regexp(S,BRACKET_POLY,'match','once');
        bracket = regexp(bracketterm,BRACKET_POLY,'match','once');
        rep = regexp(bracket,POLY,'match','once');
        S = regexprep(S,BRACKET_POLY,rep,'once');
    end
    
end