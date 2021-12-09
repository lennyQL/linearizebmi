%% Pole place probklem

% Plant
n = 4;
m = 2;
A = rand(n,n);
B = rand(n,m);
C = rand(m,n);

% Decision variables
P = sdpvar(n,n);
K = sdpvar(m,m);
% Feasible solutions
P0 = rand(n,n);
K0 = rand(m,m);
% Decomposition matrix (satisfy [He(G)>0])
G = eye(m,m);

% BMI as a string
Fstr = "P*(A+B*K*C)+(A+B*K*C)'*P";

% convert BMI into dilated LMI
[LMI,Lstr] = linearizebmi(Fstr,{'P','K'},{'P0','K0'},'G')
constraints = [LMI <= 0]