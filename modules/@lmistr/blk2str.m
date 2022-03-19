function [str, comp] = blk2str(LMIstr)
%BLK2STR 

if LMIstr.n == 5
    str = [LMIstr.Q   LMIstr.LXN;...
           LMIstr.GYR LMIstr.G];
elseif LMIstr.n == 10
    str = [LMIstr.Q   LMIstr.LXN LMIstr.O1;...
           LMIstr.O2  LMIstr.G1  LMIstr.G2;...
           LMIstr.HYR LMIstr.O3  LMIstr.H];
end
comp = LMIstr.comp;

end

