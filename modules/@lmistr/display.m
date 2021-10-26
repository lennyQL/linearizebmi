function display(LMIstr)
%DISP overloaded

[str,comp] = blk2str(LMIstr);
fprintf(['\n' '拡大LMIの文字列: ' 'He( LMI )' '\n']);
fprintf(['LMI:' '\n']);
disp(str);
fprintf(['\b' 'comp:\n    "' comp '"\n\n']);


end

