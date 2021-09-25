function [L,B] = sepbilinear(C,vlist)
% 線形項と双線形項の分離
%

smatrix = C;
Xstr = vlist{1};
Ystr = vlist{2};

% 行列
linearmatrix = smatrix; % 定数項，1次項
binearmatrix = smatrix; % 双線形項

% 一時的リスト
linearlist = {}; % 定数項，1次項
binearlist = {}; % 双線形項


for col=1:size(smatrix,1)
    % 各行ベクトル
    for row=1:size(smatrix,2)
        % 各行列要素
        % disp(col+" "+row)
        termlist = smatrix{col,row};
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

L = linearmatrix;
B = binearmatrix;

end

