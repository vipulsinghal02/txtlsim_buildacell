function medianval = curvewiseMedian(dataMat)
	% Compute the median of a set of curves given in dataMat
	% (of dimensions #points x #curves) using the max value of
	% each curve to order the curves. If there is an even number
	% of curves, then the larger curve is used. 


	[mx, I]  = sort(max(dataMat));
	if mod(numel(I), 2)==0 
		ix = I(numel(I)/2+1);
	else
		ix = I((numel(I)+1)/2);
	end
	medianval = dataMat(:,ix);
end


