
%GREEDYSCP Greedy SCP algorithm.
%	[SolC,SolL] = GREEDYSCP(C, L) if C is an array, creates a cell array SolC that is a solution of Set Cover Problem defined by C, where C{i} = S_i, an input set made by some of the elements we want to cover; SolC is made by the cells of C selected by the algorithm. The elements that we want to cover are indicates by numbers from 1 to n, where n is the number of elements we want to cover; therefore, C{i} is a vector of integers between 1 and n.
%
%	If C is a logical or numerical array of n rows, where C(j,i) > 0 iff element j is contained in set S_i, the output SolC will be a logical array made by the column of log(C) corresponding to the solution
%	
%	If a vector L of integer labels of the elements of C is provided, SolL contains the labels corresponding to SolC. Otherwise SolL contains the positions of elements of SolC in C. SolC and SolL elements are sorted in ascending order of SolL.
%
%	This is an implementation of the well-known greedy algorithm (Chvátal, 1979), with two small modifications:
%	* In case of more than one possible choice at one step, the biggest set is chosen.
%	* Once the solution is found, we check the selected sets to find a better cover solution, removing a set if is a subset of the union of the other set.
%	
%	If you use this code, please cite:
%	F. Gori, G. Folino, M.S.M. Jetten, E. Marchiori
%	"MTR: Taxonomic annotation of short metagenomic reads using clustering at multiple taxonomic ranks", Bioinformatics 2010.
%	doi = 10.1093/bioinformatics/btq649 
%
%*=*=*=*=*
%* Input adjustment
% Checking first input

function [solution_setsCell, solution_setsLabelsV] = greedyscp(setsCell, setsLabelsV)

cellL = iscell(setsCell) ;
if cellL
	if any(cellfun(@isempty, setsCell))
		disp('Hata - Giriþ dizisinde bazý boþ hücreler var')
		return;
	end
	
	setsCell = cellfun(@(x) x(:)', setsCell, 'Çýktý', false) ;

	setsN = numel(setsCell) ;
	if size(setsCell) ~= [setsN 1]
		setsCell = setsCell' ;
	end
	OriginalElementsV = unique([setsCell{:}]) ;
	elementsN = numel(OriginalElementsV) ;
	if (min(OriginalElementsV) == 1) && (max(OriginalElementsV) ==  elementsN)
		
		setsCardinalitiesV = cellfun(@numel, setsCell) ;
	else
		disp('hata ')
		return;
	end
else
	A = setsCell ;
	if isnumeric(A)
		disp('Uyarý: girdi sayýsal bir dizidir; onu mantýksal bir diziye dönüþtürmeye çalýþýyorum ')
		A = logical(A) ;
	end	
	if islogical(A)
		if ~all(sum(A, 2))
			disp('Hata - Giriþ dizisinde bazý boþ satýrlar var! ') 
			return;
		end
		[elementsN, setsN] = size(A) ;
		setsCardinalitiesV = sum(A, 1) ;
		if ~all(setsCardinalitiesV)
			disp('Hata - Girdi dizisinde bazý boþ sütunlar var!') 
			return;
		end
	else
		disp('Hata - Girdi dizisi "hücre", "mantýksal" veya "sayýsal" türünde deðil!') 
		return;
	end
end

if nargin == 1
	disp('Uyarý - Sütun etiketleri oluþturma'); 
	setsLabelsV = (1 : setsN)' ;
else 
	if size(setsLabelsV) ~= [setsN 1]
		setsLabelsV = setsLabelsV' ;
	end
end

if cellL
	[solution_setsCell, solution_setsLabelsV] = cellSCP(setsCell, setsLabelsV, setsCardinalitiesV, OriginalElementsV, elementsN, setsN) ;
else
	[solution_setsCell, solution_setsLabelsV] = matrixSCP(A, setsLabelsV, setsCardinalitiesV) ;
end

function [solution_setsCell, solution_setsLabelsV] = cellSCP(setsCell, setsLabelsV, setsCardinalitiesV, OriginalElementsV, elementsN, setsN)

numCoverVector = zeros(1, elementsN) ;

posCovSetVector = zeros(1, elementsN) ;
for iSet = 1 : setsN 
    elementsCovLogical = ismember(OriginalElementsV, setsCell{iSet}) ;
    numCoverVector = numCoverVector + elementsCovLogical ;
	
    posCovSetVector(elementsCovLogical) = iSet ;
end

UnCovLogical = numCoverVector == 1 ;
if any(UnCovLogical)
    toBeSelectVector = unique(posCovSetVector(UnCovLogical)) ;
    display([num2str(sum(UnCovLogical)) 'öðeler benzersiz bir þekilde kaplý ']);
    display([num2str(numel(toBeSelectVector)) ' setler seçilecek']);
    toBeSelectVector = unique(posCovSetVector(UnCovLogical)) ;
    solution_setsCell = setsCell(toBeSelectVector) ;
    solution_setsLabelsV = setsLabelsV(toBeSelectVector) ;  
    % güncelee
    setsCell(toBeSelectVector) = [ ] ;
    setsLabelsV(toBeSelectVector) = [ ] ;
    setsCardinalitiesV(toBeSelectVector) = [ ] ;
    %
	elementsCoveredV = unique([solution_setsCell{:}]) ;
	display(['Elements Covered: ' num2str(numel(elementsCoveredV))])
    % Equivalent: remainingElementsVector = setdiff(OriginalElementsV, elementsCoveredV) ;
    remainingElementsVector = OriginalElementsV(~ismember(OriginalElementsV, elementsCoveredV)) ;
end

%%%Açgözlü algoritmayý çalýþtýrýn
if isempty(remainingElementsVector)
	disp(' Algoritmaya gerek yok')
else
	disp(' ')
	
	[setsCardinalitiesV, sortedIndexVector] = sort(setsCardinalitiesV, 'aþaðý yuvarla') ;
	setsLabelsV = setsLabelsV(sortedIndexVector) ;
	setsCell = setsCell(sortedIndexVector) ;
    disp('Greedy algoritmasý: baþla')
	countN = 0 ;
    while ~isempty(remainingElementsVector) 

        thresholdN = sum(ismember(setsCell{1}, remainingElementsVector)) ;
        
        indexFocusedSetsLogical = setsCardinalitiesV >= thresholdN ;
        focusedSetsCell = setsCell(indexFocusedSetsLogical) ;
        focusedLabelsVector = setsLabelsV(indexFocusedSetsLogical) ;
        focusedN = sum(indexFocusedSetsLogical) ;
		