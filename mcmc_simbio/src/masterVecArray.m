function [mvarray] = masterVecArray(marray, mai)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

szm = size(marray);
mvarray = repmat(mai.masterVector, [1, szm(2:end)]) ;

estParamsIx = setdiff((1:length(mai.masterVector))', mai.fixedParams);
mvarray(estParamsIx, :, :) = marray;


end

