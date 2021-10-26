function [str, comp] = blk2str(LMIstr)
%BLK2STR 

str = [LMIstr.Q LMIstr.LXN; LMIstr.GYR LMIstr.G];
comp = LMIstr.comp;

end

