

% To run GWO: [Best_score,Best_pos,GWO_cg_curve]=GWO(SearchAgents_no,Max_iteration,lb,ub,dim,fobj)
%__________________________________________
clc;
clear;
close all;

data = load('mydata');
X = data.X;
k = 3;

CostFunction=@(m) ClusteringCost(m, X);     % Cost Function

VarSize=[k size(X,2)];  %  k kadar yani 3 e 2 lik matrix oluþturarak random m deðrini ayarlýyor 

nVar=prod(VarSize);     % Number of Decision Variables

VarMin= repmat(min(X),k,1);      % küme merkezi rasgele sýnýrlarýný atýyor 
VarMax= repmat(max(X),k,1);      % Upper Bound of Variables
VarMin=min(VarMin(:));
VarMax=max(VarMax(:));

SearchAgents_no=30; % Number of search agents

Max_iteration=100; % Maximum numbef of iterations

% Load details of the selected benchmark function
%[lb,ub,dim,fobj]=Get_Functions_details(Function_name);

[Best_score,Best_pos,GWO_cg_curve]=GWO(SearchAgents_no,Max_iteration,VarMin,VarMax,nVar,CostFunction);



%_-------------------bestSol------------------------------?
empty_individual.Position=[];
empty_individual.Cost=[];
empty_individual.Out=[];
pop=repmat(empty_individual,SearchAgents_no,1);

for i=1:SearchAgents_no
    
    % Initialize Position
    pop(i).Position=unifrnd(VarMin,VarMax,VarSize);
    
    % Evaluation
    [pop(i).Cost, pop(i).Out]=CostFunction(pop(i).Position);
    
end
% Sort Population
Costs=[pop.Cost];
[Costs, SortOrder]=sort(Costs);
pop=pop(SortOrder);
% Store Best Solution
BestSol=pop(1);
%-------------------------------------------------?



