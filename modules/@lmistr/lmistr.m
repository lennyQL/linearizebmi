function LMIstr = lmistr(Q,LXN,GYR,G,comp)
%LMISTR 
% linearizebmiの戻り値2LMIstrのオブジェクトクラス

LMIstr.Q = Q;
LMIstr.LXN = LXN;
LMIstr.GYR = GYR;
LMIstr.G = G;
LMIstr.comp = comp;

LMIstr = class(LMIstr, 'lmistr');


end

