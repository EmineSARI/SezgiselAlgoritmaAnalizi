

Choices = {'Genetic Algorithm (GA)', 'Particle Swarm Optimization (PSO)', 'Grey Wolf Algorithm(GWO)'};

ANSWER = questdlg('Select the algorithm to perform clustering.', ...
                  'Evolutionary Clutsreing', ...
                  Choices{1}, Choices{2}, Choices{3}, ...
                  Choices{1});

if strcmpi(ANSWER, Choices{1})
    ga;
    
      disM1=pdist(X);
distM=squareform(disM1);
ind=pop(1).Out.ind;
disp(sprintf('Dunns index for GA %d', dunns(k,distM,ind)));
    
    return;
end

if strcmpi(ANSWER, Choices{2})
    pso;   
    
    
    disM1=pdist(X);
distM=squareform(disM1);
ind=particle(1).Best.Out.ind;%k�meleme i�lemini yat�ktan sonra her bir veri setinin hagi k�mede olddu�unu numarland�rd���nda t�m veri setini bu de�i�kene kaydediyoruz
disp(sprintf('Dunns index for PSO %d', dunns(k,distM,ind)));%veri setinin her bir verisini k�me merkezine olan yak�nl���na g�re numara atand� �rn 5 k�me merkezi varsa 5 tane numara atarak hangi k�meye ait oldu�unubulan veriyi dunn fonk ile en iyi sonucu bulamya g�dneriyoruz
    
    
    return;
end

if strcmpi(ANSWER, Choices{3})
    mainGWO;
    
disM1=pdist(X);
distM=squareform(disM1);
ind=BestSol.Out.ind  ;
disp(sprintf('Dunns index for GWO %d', dunns(k,distM,ind)));

    return;
end;
