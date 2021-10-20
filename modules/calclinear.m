function eval = calclinear(termlist)
% 線形項の各cellの要素を計算する関数

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
    %                     if regexp(var,'(?<!\D+)\d+')
    %                         % 数値の場合
    %                         qeval = qeval * str2double(var);
    %                     else
    %                         % 変数名の場合
    %                         qeval = qeval * evalin('base', var);
    %                     end
                qeval = qeval * evalin('base', var);
            end
        end
        eval = eval + qeval;
    end

    % 要素に項がない場合，ゼロ行列を入れる
    if isempty(termlist)
        eval = zeros(colsize(col),rowsize(row));
    end

end

