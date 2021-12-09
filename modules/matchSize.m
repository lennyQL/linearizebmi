%% Useful Functions

function [M1,M2] = matchSize(L1,L2)
% 2つのリストのサイズを合わせる
% Nanを追加することで.
% グラフ出力用のデータ整形関数

sizeL1 = length(L1);
sizeL2 = length(L2);

M1 = L1;
M2 = L2;

sizediff = sizeL1 - sizeL2;

if sizediff > 0
    for i=1:abs(sizediff)
        M2 = [M2, NaN];
    end
elseif sizediff < 0
    for i=1:abs(sizediff)
        M1 = [M1, NaN];
    end
end


end