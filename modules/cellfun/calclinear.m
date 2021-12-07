function celleval = calclinear(smatrix,colsize,rowsize)
% 線形項の各cellの要素を計算する関数

celleval = smatrix;

for col=1:size(smatrix,1)
    % 各行ベクトル
    for row=1:size(smatrix,2)
        % 各行列要素
        % disp(col+" "+row)
        termlist = smatrix{col,row};
        eval = 0;
        for i=1:size(termlist,1)
    %         disp("----------")
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
            eval = eval + qeval;
        end

        % 要素に項がない場合，ゼロ行列を入れる
        if isempty(termlist)
            eval = zeros(colsize(col),rowsize(row));
        end
        
        celleval(col,row) = {eval};
        
    end
end

end

