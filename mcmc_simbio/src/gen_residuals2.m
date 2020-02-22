function res = gen_residuals2(logp, exportedMdlObj, data_array, timevec, ...
	dosedInitVals, measuredSpecies)
% NOT COMPLETE YET
% this is different from gen_residuals in that the data_array here is 
% #timepoints X #doses X #measured species. 
%
% dosedInitVals needs to be a nDoseCombinations x nSpeciesToDose matrix. 
% the concentration of this matrix needs some work from the dosedArray struct. 
%
% logp should be a column vector. in exportedMdlObj, we already have the things 
% that can be varied to be just the parameters to be estimated and the species that can be dosed, in that order. 
% measuresSpecies is a cell array of strings. 
% (unlike the struct in objFun_multiDoseSpecies.m)
% 
meanVals = mean(mean(data_array,1),2);

% meanVals = mean(reshape(data_array, [size(data_array, 1) size(data_array, 3)]),1);
wt = sum(meanVals)./meanVals; %hight mean = lower wt 
%!TODO is the mean the right statistic, or is the median or max a better statistic?
relWt = wt/sum(wt);
CONC_temp = zeros(length(timevec)*size(dosedInitVals,1), length(measuredSpecies));
for i = 1:size(dosedInitVals,1)

	sd = simulate(exportedMdlObj, [exp(logp); dosedInitVals(i,:)']);
	sd = resample(sd, timevec);

		spSD = selectbyname(sd, measuredSpecies); 
		% !TODO :test if i can just feed in the cell of strings. 
		CONC_temp((i-1)*length(timevec) + 1:(i-1)*length(timevec) +length(timevec), :) = spSD.Data; % and then just do this. 
end
ExpData = data_array;
if isequal(size(CONC_temp), size(ExpData))
			residuals = repmat(relWt, size(CONC_temp,1),1).*(CONC_temp - ExpData); 
			res = residuals(:);
            % scale the residuals by relative weight. THis emphasizes
            % the species which are low conc,
            % and deemphasizes species which are high in magnitude
            % to make them all equally important
else
    error('txtl_toolbox:objFun_multiDoseSpecies:incompatibleArrays',...
        'The simulation and experimental data need to be the same sizes.')
end

% I am not sure how to get away from a for loop. 
% !TODO: Ask AS, RMM, or Sam about code vectorization. 
% in exportedMdlObjArray, we have already set the dosing. 
% All that remains to be set is the parameters. 




% if i can somehow have a construct that only takes params as inputs, 
% and then just runs the models and subtracts from the data, instead of 
% doing the fol loop after the simulations. 
% I dont think I can do that. So lets not try. zzz

% can try interp1 or sd.resample to see which is faster. 
% I think deciding which strategy is faster will require some testing. 

% for now lets just use the naive sd.resample method. 
% If its abysmally slow, can try to modify this code to have interp1, 
% with the measured species indices predetermined. 

