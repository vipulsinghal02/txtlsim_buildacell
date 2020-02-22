function llike = gen_residuals_4(logp, em, data_array, tv, ...
    dv, ms, logresvec, stdev)
%{
% This code is an intermediary between existing code (gen_residuals_3) and
% the future code that will be the final version. Basically here I
% prototype the capability for the computation of the log likelihood a bit.
% In particular I try to reduce the size of the matrices that must be kept
in memory by computing the log likelihood in parts.
%
%
% - logp is a vector of log transformed parameter (and species) values.
%
% - em is an exported simbiology model object
%
% - da (data array) is a matlab array of numbers with dimensions:
% dim 1: has the length of tv
% dim 2: species (length is # of measured species)
% dim 3: replicates (#replicates)
% dim 4: dosing / ICs
%
%
% - tv is a time vector, just a vector of timepoints in seconds
%
% - dv is a matrix of dose vals of size  # species to
% dose x # dose combinations. We do not need to specify the names of the
% species to dose because
% that gets done when the exported model object gets made.
%
% - ms is a cell array of subcells of strings. so for example, we have
% {{'species a'}, {'species b', 'species c'}} then the first output of the
% model is the trajectory of species a, and the second output is the sum of
% the trajectores of species b and c. These two outputs will respectively
% correspond to the first column and the second column of the data array
% (for a given dose). Note that the strings 'species x' must correspond to
% species in the model object.
%
% There is also a combined optimization mode that I will attempt to build here
% soon The idea is:
% It is a way to estimate shared parameters across models when I want
% to use different data to estimate sets of parameters that overlap in
% different ways across the data.
%
% Vipul Singhal, CIT 2017
    %}
    
    meanVals = mean(mean(mean(data_array, 1), 3), 4);
    % a 1 by # measured variables array
    wt = sum(meanVals)./meanVals; %hight mean = lower wt
    relWt = wt/sum(wt); % note that relWt is a row vector.
    
    CONC_temp = zeros(length(tv), length(ms), 1, size(dv,2));
    llike = 0;
    for i = 1:size(dv,2)
        sd = simulate(em, [exp(logp); dv(:,i)]);
        sd = resample(sd, tv);
        for j = 1:length(ms)
            % COMPUTE THE SIMULATED TRAJECTORY
            % each set of measures species to sum - can remove the loop
            % if each species is individual.. in the main version have a
            % different mode
            measuredspecies = ms{j};
            spSD = selectbyname(sd, measuredspecies);
            summed_trajectories = sum(spSD.Data, 2);
            CONC_temp(:, j,1, i) = summed_trajectories;
            
            % COMPUTE THE RESIDUAL - can use repmat here because the
            % matrices are probably not big enough to slow the code down.
            % On the other hand, the for loop to do the replicates might
            % actually slow the code down.
            relWt_tiled = repmat(relWt(1,j), size(CONC_temp,1), 1,  size(data_array,3));
            replicatedsimdata = repmat(CONC_temp(:, j,1, i), [1, 1, size(data_array, 3)]);
            residuals = relWt_tiled.*(replicatedsimdata - data_array(:, j, :, i));
            res = residuals(:);
            llike = llike + sum(logresvec(res, stdev));
        end
    end
    
end




