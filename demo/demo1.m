%% 極配置
n = 4;
m = 2;
A = rand(n,n);
B = rand(n,m);
C = rand(m,n);

% 決定変数
P = sdpvar(n,n);
K = sdpvar(m,m);
% 暫定解
P0 = rand(n,n);
K0 = rand(m,m);
% 分割行列
G = eye(m,m);

Fstr = "P*(A+B*K*C)+(A+B*K*C)'*P";

[LMI,Lstr] = linearizebmi(Fstr,{'P','K'},{'P0','K0'},'G')
constraints = [LMI <= 0]