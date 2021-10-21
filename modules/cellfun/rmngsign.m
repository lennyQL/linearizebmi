function smatrix = rmngsign(termlist)
% 項のマイナス(-)符号を削減

smatrix = {};

for i=1:size(termlist,1)
    % disp("----------")
    term = termlist{i,1};
    varlist = {};
    signcount = 0; % -符号のカウンタ
    
    % 符号探索
    for j=1:size(term,2)
        var = term{1,j};
        if var == "-"
            signcount = signcount + 1;
        else
            varlist = updateList(varlist,var);
        end
    end
    
    % -が奇数個なら，項の先頭に-を付ける
    if mod(signcount,2) ~= 0
        varlist = cat(2,{["-"]},varlist);
    end
    
    % 項のリストを結合(縦のcell配列として)
    smatrix = updateList(smatrix,varlist,1);
    
end


end

