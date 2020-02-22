function llike = gen_residuals_v2(logp, estParamIx, fixedMasterVec, data_array,...
                            timeVec, mi, logresvec, stdev)
% This code is for computing the log likelihood with parameters spread out over
% multiple models (network topologies) - geometries. 
%{ OLD DOCUMENTATION from gen_residuals_4
% This code is an intermediary between existing code (gen_residuals_3) and
% the future code that will be the final version. Basically here I
% prototype the capability for the computation of the log likelihood a bit.
% In particular I try to reduce the size of the matrices that must be kept
% in memory by computing the log likelihood in parts.
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
% - mspecies is a cell array of subcells of strings. so for example, we have
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
    
    % the unpacking happens in steps. (quite similar to integrableLHS_v2)
    
fixedMasterVec(estParamIx) = logp;
fullMasterVec = fixedMasterVec;






    llike = 0;
    
    for kk = 1:length(mi) % for each topo
        if isfield(mi(kk), 'experimentWeighting')
            if ~isempty(mi(kk).experimentWeighting)
                % The relative importance of this topology 
                topoWeight = mi(kk).experimentWeighting; % a scalar number
            else 
                topoWeight = 1;
            end
        else
            topoWeight = 1;
        end
        
        
        pmaps = mi(kk).paramMaps;
        
        % ds = struct('names', {mi(kk).dosednames},...
        %          'dosematrix', mi(kk).dosedvals);
        
        dv = mi(kk).dosedVals; % can doses be reordered in the export process? 
%         This is very important to check. 
        if isfield(mi(kk), 'doseWeighting')
            % the dose weighting muse have size 1 by number of dose
            % combintations. 
            if isequal(size(mi(kk).doseWeighting), [1 size(dv,2)])  
                
                % The ralative importance of this topology is given by
                doseWeight = mi(kk).doseWeighting;
            else
                % otherwise, all the doses are equally weighted. 
                doseWeight = ones(1,size(dv,2));
            end
        else
            doseWeight = ones(1,size(dv,2));
        end
        
        em = mi(kk).emo;
        mspecies = mi(kk).measuredSpecies;
        % for each geom
        for hh = 1:size(pmaps, 2)
            
            % pmaps is in the order defined by the unordered list of names
            % in each model info. we need to reorder these indices. 
            pIX_tg = pmaps(mi(kk).orderingIx, hh); 
            % THIS REORDERING STEP IS VERY IMPORTANT. 
            
            pvec_tg = fullMasterVec(pIX_tg);
            da = data_array{kk}(:, :, :, :, hh);
            tv = timeVec{kk};

            meanVals = mean(mean(mean(da, 1), 3), 4);
            % a 1 by # measured variables array
            wt = sum(meanVals)./meanVals; %hight mean = lower wt
            relWt = wt/sum(wt); % note that relWt is a row vector.
            
            CONC_temp = zeros(length(tv), length(mspecies), 1, size(dv,2));
           
            for ii = 1:size(dv,2)
                
                % pvec_tg needs to be in the ordered state, ie,
                % mi(kk).namesOrd. 

                %try
                    sd = simulate(em, [exp(pvec_tg); dv(:,ii)]);
                %catch ME
                %    disp(ME.identifier);
                %end
                
                sd = resample(sd, tv);
                for jj = 1:length(mspecies)
                    % COMPUTE THE SIMULATED TRAJECTORY
                    % each set of measures species to sum - can remove the loop
                    % if each species is individual.. in the main version have a
                    % different mode
                    measuredspecies = mspecies{jj};
                    spSD = selectbyname(sd, measuredspecies);
                    summed_trajectories = sum(spSD.Data, 2);
                    CONC_temp(:, jj,1, ii) = summed_trajectories;
                    
                    % COMPUTE THE RESIDUAL - can use repmat here because the
                    % matrices are probably not big enough to slow the code down.
                    % On the other hand, the for loop to do the replicates might
                    % actually slow the code down.
                    relWt_tiled = repmat(relWt(1,jj), size(CONC_temp,1), 1,  size(da,3));
                    replicatedsimdata = repmat(CONC_temp(:, jj,1, ii), [1, 1, size(da, 3)]);
                    residuals = relWt_tiled.*(replicatedsimdata - da(:, jj, :, ii));
                    
                    % multiply the residuals with the topology's relative
                    % importance, and the dose's relative importance. 
                    res = topoWeight*doseWeight(ii)*residuals(:);
                    llike = llike + sum(logresvec(res, stdev));
                end
            end
        end
    end
    
        
end




