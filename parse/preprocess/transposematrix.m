
%% 行列の転置([...;...])'を展開する関数，文字列のままで
function S = transposematrix(S)

%     regdeclare();
    global MAT_TRANSPOSE
    global VEC
    global MAT_COL
    global POLY_BRACKET
    trans = regexp(S,MAT_TRANSPOSE,'match','once');
    matrix = regexp(trans,VEC,'match','once');
    col = regexp(matrix,MAT_COL,'match');
    item = regexp(col,POLY_BRACKET,'match');

    for i=1:length(item)
        if i == 1
            itemlist = item{1,i};
        else
            itemlist = cat(1,itemlist,item{1,i});
        end 
    end
    % itemlist
    % itemlist = itemlist'
    itemlistT = "("+itemlist'+")'";

    smat = "[";
    for i=1:size(itemlistT,1)
        col = strjoin(itemlistT(i,:)," ");
        smat = smat+col;
        if i ~= size(itemlistT,1)
            smat = smat+";";
        else
            smat = smat+"]";
        end 
    end
    
    S = regexprep(S,MAT_TRANSPOSE,smat,'once');
    S = char(S);

end
