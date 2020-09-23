

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
ind=particle(1).Best.Out.ind;%kümeleme iþlemini yatýktan sonra her bir veri setinin hagi kümede oldduðunu numarlandýrdýðýnda tüm veri setini bu deðiþkene kaydediyoruz
disp(sprintf('Dunns index for PSO %d', dunns(k,distM,ind)));%veri setinin her bir verisini küme merkezine olan yakýnlýðýna göre numara atandý örn 5 küme merkezi varsa 5 tane numara atarak hangi kümeye ait olduðunubulan veriyi dunn fonk ile en iyi sonucu bulamya gödneriyoruz
    
    
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
