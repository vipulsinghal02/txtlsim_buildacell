function outputslice = mcmc_cut(marray, pID, pRanges, varargin)
% Cut a lower dimensional slice out of a high dimensional
% parameter distribution. 
% 
% - marray: This function takes an array of numbers of size 
% nPoints x nParam, and returns an array of size 
% nPointsFewer x nParam, where nPointsFewer is the number of 
% points satisfying the slicing condition described below. 
% marray can also be an array of size nParam x nWalkers x nsteps, 
% in which case the output is of the same dimension, or of dimension
% nParam x nWalkersSmaller, where nWalkersSmaller is the total 
% number of points where the slicing condition holds. Note that the walker
% positions with respect to the steps is not preserved.
%
% - pID is a vector of parameter indices ranging from 1 to 
% nParam. 
% 
% - pRanges is a 2 x length(pID) array of upper and lower bounds 
% (first and second row resp) for parameter values corresponding 
% to the parameters indexed by the array pID. 
%
% SLICING CONDITION:
% This function returns the subset of points (rows) in marray
% for which the condition all(pranges(1, :) >= point(pID)) and 
% all(pranges(2, :) <= point(pID)), where the vector 'point'
% is a row of the array of parameter points marray. (The Output
% 'outputslice' is a matrix of all such points.)
% 
% Optional Input name-value pair
% 'marginalize' (default: false). If set to true, we remove 
% the parameters with indices in the index vector pID from 
% the returned array, and and in this case the output array 
% OUTPUTSLICE will have dimensions 
% nPointsFewer x (nParam - length(pID)). 
% 
p = inputParser;
p.addParameter('marginalize', false, @islogical)
p.parse(varargin{:});
p = p.Results;
convertTo3D = false;
if ndims(marray) == 3
    convertTo3D = true;
    szm = size(marray);
    marray = marray(:, :)';
end

fun1 = @(rowvec) all(rowvec(pID)<pRanges(1, :)) &...
				 all(pRanges(2, :)<rowvec(pID));

indices = zeros(size(marray, 1), 1);
for i = 1:size(marray, 1)
	indices(i) = fun1(marray(i, :));
end

if p.marginalize
	ixtoreturn = setdiff(1:size(marray, 2), pID);
	outputslice = marray(indices, ixtoreturn);
else
    if sum(indices) == 0
        outputslice = [];
        warning('parameter ranges too narrow to capture any points.');
    else
        
	outputslice = marray(logical(indices), :);
    end
    
end

    if convertTo3D
        nsteps = floor(size(outputslice, 1)/szm(2));
        outputsliceold = outputslice';
        if ~(nsteps == 0)
            npoints = szm(2)*nsteps;
            outputslice = reshape(outputsliceold(:,(end-npoints+1):end),...
                [szm(1), szm(2), nsteps]);
        end
    end

end

