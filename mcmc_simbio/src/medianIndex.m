function [ix, medianvals] = medianIndex(inputarray, dim)
	% compute the index of the median element of an
	% array along a specified dimension. If there is an
	% even number of elements, then pick the bigger of the middle 
	% two elements. 
	% ix is an array of the same size as 'inputarray', but with a width 
	% of 1 on dimension number dim. The elements of ix give the element 
	% along 'dim' in the array 'inputarray' that is either the median 
	% element along that dimension, or if there are an even number of 
	% element in that dimension, then is the closest overestimate of 
	% the median. The output medianvals is an array of the values in 
	% inputarray that ix points to. 
	% 
	% (c) Vipul Singhal
    
    % EXAMPLE: inputarray is, for eg, 1-2-1-4 ---> 1-4-1-2 when dim = 3. 
    % time - ms - replicates - doses ---> replicates - doses - time - ms
	shiftedarray = shiftdim(inputarray, dim-1);
    
   
    % if all the leading dimensions to be shifted are singletons, then they
    % will be lost when shiftdim works: ie, if dim = 3, so the shifting is
    % to be by 2 dimensions, then the result of shifting is: 
    % 1-2-3-4 ---> 3-4-1-2 (stays 4D array!)
    % 1-1-3-4 ---> 3-4 (becomes 2D array!)
    % 
    % Also, note that 
    %     rr = rand(2, 1, 3, 4);
    %     size(shiftdim(rr, 2))
    % 
    %     ans =
    % 
    %          3     4     2
    % 
    % ie, basically matlab will just not report trailing singletons, and
    % therein lies the problem. 
    % 
    % 
    % so later when we do the shift back, have to do it as follows: 
    % 
    % First pad on the left with singletons equal to the number removed:
    % 
    % [1-2-3-4 --->] 3-4-1-2 --(pad with 0)--> 3-4-1-2, then rotate by 
    % ndims - ((ndims - dim +1) - 0) == 2, ---> 1-2-3-4. [ORIGINAL]
    % 
    % 
    %
    % [1-1-3-4 --->] 3-4 --(pad with 2 singletons)--> 1-1-3-4, then rotate
    % by ndims - ((ndims - dim +1) - 2) == 4, ---> 1-1-3-4  [ORIGINAL]
    %
    %
    %
    % [2-1-3-4 --->] 3-4-2 --(pad with one singleton)--> 1-3-4-2 
    % then rotate by ndims -((ndims - dim +1) - 1) == 3, ---> 2-1-3-4 [ORIGINAL]
    % 
    % Then do the rotation in the same direction as the original, but 
    % ndims - dim +1
    % 

    %{ OLD:
    % This can be undone by first rotating the 
    % remaining dimensions back, and then adding singletons on the right. 
    %     if ndims(shiftedarray) < ndims(inputarray)
    %         dimsToRightPad = ndims(inputarray) - ndims(shiftedarray);
    %     end
       % otherwise: 
	% undo using shiftdim(shiftedarray, ndims(array) - dim +1)
    % actually instead of circularly shifting, can just do :
    %     shiftdim(I5, -dim+1) to undo it!!!! <--- NOPE. 
    % From the documentation: "When N is negative, shiftdim
    % shifts the dimensions to the right and pads with singletons."
    %}
    
    % EXAMPLE (cont): srted is 1-4-1-2, and in this case the same as 
    % shiftedarray, since dimension 1 has length 1, so there is nothing to
    % sort. I is a array of ones of the saem size. 
	[srted, I] = sort(shiftedarray, 1);

	%index of the element in srted that is the closest 
	% overestimate of the median 
	% also indices the first dimension of I
    
	II = ceil((size(I, 1)+1)/2); 

	% Reshape I into a matrix
	III = I(:,:);
	srted3 = srted(:,:);
	% index using II (ie, find the closest-to-middle element)
	I4 = III(II,:);
	srted4 = srted3(II, :);
	% Now reshape back using the original dimensions

	szmat = size(I);
	I5 = reshape(I4, [1, szmat(2:end)]);
	srted5 = reshape(srted4, [1, szmat(2:end)]);
	% shift back to the original dimensions
    % there are 2 cases here: 
    % if the original count of the dimension sizes was 
    % 1-1-nRep-nDoses, then the shiftdim above resulted in a 2D matrix
    % 
    
    singletonsToPad = ndims(inputarray) - ndims(shiftedarray);
    padded_srted5 = shiftdim(srted5, -singletonsToPad);
    padded_I5 = shiftdim(I5, -singletonsToPad);
    
    % dimsToRotateBy = ndims -((ndims - dim +1) - singletonsToPad) 
    % = dim -1 + singletonsToPad
    dimsToRotateBy = dim - 1 + singletonsToPad;
	ix = shiftdim(padded_I5,dimsToRotateBy);
	medianvals = shiftdim(padded_srted5,dimsToRotateBy);

	


% coding notes:
	% first bring the dimension along which the median is to
	% be found to the first dimension.
	%%%%	
	% dim is 3
	% 4 total # of dims. 
	% 1 2 3 5 (sumovertime, ms, rep, dose)
	% > shiftdim by dim - 1
	% shifted has size: 3 5 1 2 (rep, dose, sumovertime, ms)
	% ie, 3 element columns, and 5 columns per pane, 1 pane per set1, 
	% 2 set1's. 
	%%%%


	% then sort along this dimension, and subsequently pick 
	% out the (possibly closest overestimate of the) median value.

	%%%%
	% sort 
	% pick out median indices

	% Ix has size 1 5 1 2 and points to elements of shifted
	% shifted(Ix) is a 1 5 1 2 array of the median elements. 
	% more than that, Ix is what we want.
	% in fact, unshifting Ix as follows: 
	% shiftdim(Ix, numel(size(array)) - dim + 1) brings Ix to 
	% 1 2 1 5 <---- Ix2
	% Now use Ix2 to index the full array back in the calling function 
	% array2(:, Ix(1,:,:,:)) maybe? I have no idea. Try it out I guess. 
	% no i dont think this works. 




	% % bring dim to be the first dimension
	% shftsrted = shiftdim(srted, dim-1); 
	% medianarray = shftsrted(ceil((size(I, dim)+1)/2), )



	%%%%







	end
