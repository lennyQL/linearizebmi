function Qchar = linear2str(Q,colsize,rowsize)
% 線形項Qの文字列を返す関数(cell配列を受け取る)

% Qの文字列
Qchar = [];

for col=1:size(Q,1)
    % 各行
    
    Qcharcol = [];
    for row=1:size(Q,2)
        % 各列
        termlist = Q{col,row};
        Qcharrow = [];
        for i=1:size(termlist,1)
            % 各要素
            term = termlist{i,1};
            qchar = [];
            for j=1:size(term,2)
                % 各項
                var = term{1,j};
                if isempty(qchar)
                    qchar = [qchar char(var)];
                elseif qchar == '-'
                    % qchar = ['(' qchar char(var) ')'];
                    qchar = [qchar char(var)];
                else
                    qchar = [qchar '*' char(var)];
                end
            end
            
            % 要素内の項を+で統合
            if isempty(Qcharrow)
                Qcharrow = [Qcharrow qchar];
            else
                Qcharrow = [Qcharrow '+' qchar];
            end
        end

        % 要素に項がない場合，ゼロ行列を入れる
        if isempty(termlist)
            var = [func2str(@zeros) '(' num2str( colsize(col) ) ',' num2str( rowsize(row) ) ')'];
            Qcharrow = [Qcharrow  var];
        end
        
        % 列の分割
        Qcharcol = [Qcharcol string(Qcharrow)];
        
    end
    % 行の分割
    Qchar = [Qchar ; Qcharcol];
end


end

