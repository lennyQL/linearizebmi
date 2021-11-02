%%% SDP solver の設定
% solversetup;

%%% 問題の定義
% 制御対象データの定義
% Complieb のデータのロード
prob='HE1';
[a,b1,b2,c1,c2,d11,d12,d21,na,nw,nu,nz,ny]=COMPleib(prob);

% 決定変数の定義
p=sdpvar(na,na,'symmetric');	% Lyapunov 行列
k=sdpvar(nu,ny,'full');		% 制御器(定数ゲイン)
g=sdpvar(1,1);			% H∞ノルム

% BMI 最適化問題の定義
Fstr = "[p*(a+b2*k*c2)+(p*(a+b2*k*c2))',p*(b1+b2*k*d21),(c1+d12*k*c2)';"   +...
        "(p*(b1+b2*k*d21))',            -g*eye(nw,nw),     (d11+d12*k*d21)';" +...
        "c1+d12*k*c2,                   d11+d12*k*d21,  -g*eye(nz)]";

% 暫定解の宣言
p0=sdpvar(na,na);
k0=sdpvar(nu,ny);
%p0=rand(size(p));
%k0=rand(size(k));

[LMIauto,Lstr]=linearizebmi(Fstr,{'p','k'},{'p0','k0'}) % 修正
LMI = [LMIauto <= 0] % 追加

% linearizebmi
%
p0init=rand(size(p));
k0init=rand(size(k));


lcmax=200;
for lc=1:lcmax
  extLMI=LMI;
  extLMI=replace(extLMI,p0,p0init);
  extLMI=replace(extLMI,k0,k0init);

  optimize(extLMI,g)

  p0init=double(p);
  k0init=double(k);
end