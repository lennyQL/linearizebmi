%% ベクトルの積[...;...]*[... ...]を展開する関数
function S = prodvec(S)
    
%     regdeclare();
    global VEC
    global PROD_VEC
    global VEC_PROD_VEC
    global POLY

    
    if regexp(S,VEC_PROD_VEC)
        vecprodterm = regexp(S,VEC_PROD_VEC,'match','once');
        left = regexp(vecprodterm,VEC,'match','once');
        right = regexp(vecprodterm,PROD_VEC,'match','once');
        lexp = regexp(left,POLY,'match');
        rexp = regexp(right,POLY,'match');
        l = "("+lexp'+")";
        r = "("+rexp+")";
        matrix = append(l,"*",r);
%         matrix'
        rep = "[";
        for i=1:size(matrix,1)
            col = strjoin(matrix(i,:)," ");
            rep = rep+col;
            if i ~= size(matrix,1)
                rep = rep+";";
            else
                rep = rep+"]";
            end 
        end
    end
    S = regexprep(S,VEC_PROD_VEC,rep,'once');
    S = char(S);


end