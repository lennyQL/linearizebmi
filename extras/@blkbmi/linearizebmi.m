function [LMI,LMIstr,gBMI,BMI] = linearizebmi(X,vlist,v0list,G)

[n,m] = size(X.blocks);

if n==m
    % Perform symmetric extension
    for i = 1:n
        for j = i+1:n
            if isempty(X.blocks{i,j})
                X.blocks{i,j} = "("+X.blocks{j,i}+")'";
            elseif isempty(X.blocks{j,i})
                X.blocks{j,i} = "("+X.blocks{i,j}+")'";
            end
        end
    end
else
    error('It must be square matrix.');
end


% cell list to a constraint as string
S = cellfun(@string,X.blocks);
C = {};
for i=1:m
    C = cat(1,C,strjoin(S(i,:)));
end
C = strjoin(C,";");
C = "["+C+"]";


% output
v0listtmp = v0list;

if nargin == 1
    LMI = S;
elseif nargin == 3
    [LMI,LMIstr,gBMI,BMI] = linearizebmi(C,vlist,v0listtmp);
elseif nargin == 4
    [LMI,LMIstr,gBMI,BMI] = linearizebmi(C,vlist,v0listtmp,G);
else
    error('nargin must be 1 or 3 or 4.')
end


