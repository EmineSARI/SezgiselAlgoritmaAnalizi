
%___________________________________________________________________%
%  Grey Wold Optimizer (GWO) source codes version 1.0               %
%                                                                   %
%  Developed in MATLAB R2011b(7.13)                                 %
%                                                                   %
%  Author and programmer: Seyedali Mirjalili                        %
%                                                                   %
%         e-Mail: ali.mirjalili@gmail.com                           %
%                 seyedali.mirjalili@griffithuni.edu.au             %
%                                                                   %
%       Homepage: http://www.alimirjalili.com                       %
%                                                                   %
%   Main paper: S. Mirjalili, S. M. Mirjalili, A. Lewis             %
%               Grey Wolf Optimizer, Advances in Engineering        %
%               Software , in press,                                %
%               DOI: 10.1016/j.advengsoft.2013.12.007               %
%                                                                   %
%___________________________________________________________________%

% Gri Kurt Optimize Edici
function [Alpha_score,Alpha_pos,Convergence_curve]=GWO(SearchAgents_no,Max_iter,lb,ub,dim,fobj)



% Alpha, beta ve delta_pos ��elerini ba�lat
Alpha_pos=zeros(1,dim);
Alpha_score=inf; %B�y�tme sorunlar� i�in bunu% -inf olarak de�i�tirin

Beta_pos=zeros(1,dim);
Beta_score=inf; %maksimizasyon problemleri i�in bunlar� de�i�tirin

Delta_pos=zeros(1,dim);
Delta_score=inf; %maximization problemleri i�in bunlar� de�i�tirin

%Arama arac�lar�n�n konumlar�n� ba�lat
Positions=initialization(SearchAgents_no,dim,ub,lb);

Convergence_curve=zeros(1,Max_iter);

l=0;% d�ng� say�c�

% Main d�ng�s�
while l<Max_iter
    
    for i=1:size(Positions,1)  
       % Arama alan�n�n s�n�rlar�n�n �tesine ge�en arama arac�lar�n� geri d�n
        Flag4ub=Positions(i,:)>ub;
        Flag4lb=Positions(i,:)<lb;
        Positions(i,:)=(Positions(i,:).*(~(Flag4ub+Flag4lb)))+ub.*Flag4ub+lb.*Flag4lb;               
        
        % Her arama arac�s� i�in nesnel i�levi hesapla
        fitness=fobj(Positions(i,:));
        
        % G�ncelle Alpha, Beta, and Delta
        if fitness<Alpha_score 
            Alpha_score=fitness; %  alpha g�ncelle
            Alpha_pos=Positions(i,:);
        end
        
        if fitness>Alpha_score && fitness<Beta_score 
            Beta_score=fitness; %  betay� g�ncelle
            Beta_pos=Positions(i,:);
        end
        
        if fitness>Alpha_score && fitness>Beta_score && fitness<Delta_score 
            Delta_score=fitness; %  deltay� g�nceller
            Delta_pos=Positions(i,:);
        end
    end
    
    
    a=2-l*((2)/Max_iter); %a do�rusal olarak fron 2'den 0'a d��er
    
    % Omegas dahil arama arac�lar�n�n konumunu g�ncelleme
    for i=1:size(Positions,1)
        for j=1:size(Positions,2)     
                       
            r1=rand(); %r1, [0,1] 'de rastgele bir say�d�r
            r2=rand(); % r2  [0,1] de rastgele bir say�d�r
            
            A1=2*a*r1-a; % denklem (3.3)
            C1=2*r2; % denklem (3.4)
            
            D_alpha=abs(C1*Alpha_pos(j)-Positions(i,j)); % deklem (3.5) -1.par�a
            X1=Alpha_pos(j)-A1*D_alpha; % denklem (3.6)-1.par�a
                       
            r1=rand();
            r2=rand();
            
            A2=2*a*r1-a; % Equation (3.3)
            C2=2*r2; % Equation (3.4)
            
            D_beta=abs(C2*Beta_pos(j)-Positions(i,j)); % Equation (3.5)-part 2
            X2=Beta_pos(j)-A2*D_beta; % Equation (3.6)-part 2       
            
            r1=rand();
            r2=rand(); 
            
            A3=2*a*r1-a; % denklem (3.3)
            C3=2*r2; % denklem (3.4)
            
            D_delta=abs(C3*Delta_pos(j)-Positions(i,j)); % Equation (3.5)-part 3
            X3=Delta_pos(j)-A3*D_delta; % Equation (3.5)-part 3             
            
            Positions(i,j)=(X1+X2+X3)/3;% Equation (3.7)
            
        end
    end
    l=l+1;    
    Convergence_curve(l)=Alpha_score;
     
end













