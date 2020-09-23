%
% Copyright (c) 2015, Yarpiz (www.yarpiz.com)
% All rights reserved. Please read the "license.txt" for license terms.
%
% Project Code: YPML101
% Project Title: Evolutionary Clustering in MATLAB
% Publisher: Yarpiz (www.yarpiz.com)
% 
% Developer: S. Mostapha Kalami Heris (Member of Yarpiz Team)
% 
% Contact Info: sm.kalami@gmail.com, info@ yarpiz.com
%

function [z, out] = ClusteringCost(m, ~)
     data = load('mydata');
     X = data.X;
     k=3;
    % Uzakl�k Matrisini Hesapla
    m=reshape(m,[k size(X,2)]);
    d = pdist2(X, m); %baslang��ta k�me merkezleri random
    
    % K�me Atama ve En Yak�n Mesafeleri Bulma
    [dmin, ind] = min(d, [], 2);
    
    % K�me ��indeki Mesafe Toplam�
    WCD = sum(dmin); % dmin lerin toplam�ndan daha optimum bir merkez olu�turuyor
    z=WCD;
    out.d=d;
    out.dmin=dmin;
    out.ind=ind;
    out.WCD=WCD;
    
end