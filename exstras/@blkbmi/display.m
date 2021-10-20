function F = display(X)
% DISPLAY Overloaded

try
    X = bmiparser(X);
    display(X);
catch
    disp('Incomplete block variable.');
end


