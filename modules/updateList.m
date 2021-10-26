%% Listの更新（追加）と初期化
function [L,V] = updateList(L,V,n)
    % ex) 
    %   LにVを追加，Vを初期化:
    %     [L,V] = updateList(L,V)
    %   LにVを追加のみ: 
    %     L = updateList(L,V)

    if nargin<3
        n = 2;
    end

    % Listの追加
    % n=1: 縦に追加
    % n=2: 横に追加
    L = cat(n,L,{V});
    
    if isa(V, "string")
        % stringの初期化
        V = "";
    elseif isa(V, "cell")
        % cellの初期化
        V = {};
    end
end


