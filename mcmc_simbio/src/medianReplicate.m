function medianCurves = medianReplicate(dataArray, ix)
	% Pick out the median replicate from data array, the replicate to be picked 
	% out specified by the index array ix. 
	% 
	% This function takes two inputs: dataArray and ix. dataArray is an array 
	% containing data that has dimensions time x measured species x replicates
	% x doses. ix is an array of dimensions 1 x nMS x 1 x nDoses OR 
	% nTimePoints x nMS x 1 x nDoses. If it is the latter, only the element 
	% corresponding to the first time point is used. 
	% I.e., for the ith measures species, and the jth dose, the index of the 
	% replicate within dataArray that is used comes from ix(1,i, 1,j), so that
	% if ix(2:end,i, 1,j) exist, they are not used in any way. The function that
	% creates the ix array is medianindex, and this function is used, for example
	% the computeDataStats function above, or in mcmc_runsim.m. 


		% data array must be time x measuredspecies x replicates x doses
  
        % Get the median curves. Can this be done via pure vector indexing? 
        % !UNIMPORTANT but could be fun to think about sometime. 
        medianCurves = zeros(size(dataArray, 1), size(dataArray, 2), ...
        	1, size(dataArray, 4));

        for i = 1:size(dataArray, 2)
        	for j = 1:size(dataArray, 4)
        		medianCurves(:,i,1,j) = dataArray(:, i, ix(1,i, 1,j), j);
        	end
        end
end