function smatrix = rmngsign(termlist)
% ฬ}CiX(-)๐ํธ

smatrix = {};

for i=1:size(termlist,1)
    % disp("----------")
    term = termlist{i,1};
    varlist = {};
    signcount = 0; % -ฬJE^
    
    % T๕
    for j=1:size(term,2)
        var = term{1,j};
        if var == "-"
            signcount = signcount + 1;
        else
            varlist = updateList(varlist,var);
        end
    end
    
    % -ช๏ยศ็Cฬๆชษ-๐tฏ้
    if mod(signcount,2) ~= 0
        varlist = cat(2,{["-"]},varlist);
    end
    
    % ฬXg๐(cฬcellz๑ฦตฤ)
    smatrix = updateList(smatrix,varlist,1);
    
end


end

