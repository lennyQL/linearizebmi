function mat = cell2mat(C)
% function as cell to matrix
% ovarlap original cell2mat
% class accept sdpvar or double
    mat = [];
    for col=1:size(C,1)
        % 各行ベクトル
        colvec = [];
        for row=1:size(C,2)
            % 各行列要素
            colvec = cat(2,colvec,C{col,row});
        end
        mat = cat(1,mat,colvec);
    end
end