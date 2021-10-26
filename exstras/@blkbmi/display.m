function display(X)
% DISPLAY Overloaded

try
    X = linearizebmi(X);
    disp(X);
catch
    disp('Incomplete block variable.');
end


