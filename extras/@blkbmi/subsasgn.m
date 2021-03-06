function X = subsasgn(X,I,Y)
%SUBASGN (overloaded)

try
    if strcmp('()',I.type)
        if any(I.subs{1} <=0)
            error('Index into matrix is negative or zero.');
        end
        i = I.subs{1};
        j = I.subs{2};
        X.blocks{i,j} = Y;
        
    else
        error('Reference type not supported');
    end
    
catch
    error(lasterr)
end

