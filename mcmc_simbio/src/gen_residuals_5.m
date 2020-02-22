function llike = gen_residuals_5(logp, em, data_array, tv, ...
    dv, ms, logresvec, stdev, parametermap)

%{
    % Actually the simplest thing I can do rigth now is just shared CSPs. 
    % this pretty much follows the function log_likelihood_sharedCSP.m
    % so 1 model topo. and two geometries (extracts)... or actually any number
    % of geometries. and a spec for which
    params are the CSPs, and which are the extract specific parameters and
    extract specific species. 
    %
    
    THIS VERSION OF THE CODE IS THE FIRST CUT AT BUILDING MODELS WHERE
    PARAMETERS ARE SHARED ACROSS GEOMETRIES AND TOPOLOGIES
% This version of the code is not the one where we have fully general paraneter sharing
across topologies and geometries (that comes later, and I suspect will be
slower). Here we simply have up to 2 topologies (calibration and test) and
2 geometries (before going to the full arbitrary sharing, we will
generalize this to arbitrary # of topos and geos, but still in the crossed
calib-corr method. 
    
    
% in the more general code: em is an array of exported model objects. Here
we just have em1 and em2 for the two topologies. 
%    
% I will start with considering the following modes: 
    x0b, x2b, x0, 
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
% dim 5: extracts (geometries)
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
    
    
    
    espIX = parametermap{1};
    cspIX = parametermap{2};
    nESP = length(espIX); % the ESP indices in the model (not in logpjoint)
    nCSP = length(cspIX); % the CSP indices in the model (not in logpjoint)
    
    
    nEnv = size(data_array, 5);
    meanVals = mean(mean(mean(mean(data_array, 1), 3), 4), 5);
    % a 1 by # measured variables array
    wt = sum(meanVals)./meanVals; %hight mean = lower wt
    relWt = wt/sum(wt); % note that relWt is a row vector.
    
    CONC_temp = zeros(length(tv), length(ms)); 
    % dont need the other dimensions! remove from gen_residuals_4 too. 
    
    
    paramvec = zeros(nESP+nCSP,1);
    cspindices = (nESP*nEnv+1):length(logp);
    logpcsp = logp(cspindices);
    paramvec(cspIX) = logpcsp;
    
    
    llike = 0;
    for envID = 1:nEnv
        % pick out the relevant parameters from the joint parameter vector
        espindices = (envID-1)*nESP + (1:nESP);
        logpesp = logp(espindices);
        paramvec(espIX) = logpesp;
        
        
        % set the vector of parameters and species that get estimated (ie
        % non dosing values to simulate with)
        
        for i = 1:size(dv,2)
            sd = simulate(em, [exp(paramvec); dv(:,i)]);
            sd = resample(sd, tv);
            for j = 1:length(ms)
                % COMPUTE THE SIMULATED TRAJECTORY
                % each set of measures species to sum - can remove the loop
                % if each species is individual.. in the main version have a
                % different mode
                measuredspecies = ms{j};
                spSD = selectbyname(sd, measuredspecies);
                summed_trajectories = sum(spSD.Data, 2);
                CONC_temp(:, j) = summed_trajectories;

                % COMPUTE THE RESIDUAL - can use repmat here because the
                % matrices are probably not big enough to slow the code down.
                % On the other hand, the for loop to do the replicates might
                % actually slow the code down.
                relWt_tiled = repmat(relWt(1,j), size(CONC_temp,1), 1,  size(data_array,3));
                replicatedsimdata = repmat(CONC_temp(:, j), [1, 1, size(data_array, 3)]);
                residuals = relWt_tiled.*(replicatedsimdata - data_array(:, j, :, i,envID));
                res = residuals(:);
                llike = llike + sum(logresvec(res, stdev));
            end
        end
    end
    
end




