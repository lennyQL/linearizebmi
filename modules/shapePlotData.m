%% Useful Functions

function varout = shapePlotData(varargin)
% 複数のリストのサイズを合わせる
% Nanを追加することで.
% グラフ出力用のデータ整形関数

varargn = length(varargin);
lenlist = [];
for i=1:varargn
    len = length(varargin{i});
    lenlist = [lenlist,len];
end
maxlen = max(lenlist);


output = [];
for i=1:varargn
    M = varargin{i};
    sizediff = maxlen - lenlist(i);
    if sizediff > 0
        for n=1:abs(sizediff)
            M = [M, NaN];
        end
    end
    output = [output; M];
end
varout = output';

end