function res = gen_residuals_3(logp, em, data_array, tv, ...
    dv, ms, varargin)
% generate residuals for MCMC. 
% - logp is a vector of log transformed parameter (and species) values. 
%
% - em is an exported simbiology model object
%
% - da (data array) is a matlab array of numbers with dimensions of size
% #timepoints x #measured variables x #ICs (aka dose combinations)
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
% There is also a combined optimization mode that will exist in the future. 
% I have not really completed it yet, so ignore that for now. The idea is: 
% It is a way to estimate shared parameters across models when I want
% to use different data to estimate sets of parameters that overlap in
% different ways across the data. 
%
% Vipul Singhal, CIT 2017

p = inputParser;
% estimation structure is a matrix that specifies which parameters in logp
% to apply to the individual exported model objects. Ie, which parameters
% in logp get estimated for each (model, data, timevec, dose, measured
% species) tuple.
p.addParameter('multiopt_params', [], @isnumeric);
p.parse(varargin{:})
p = p.Results;
multiparam = p.multiopt_params; 
summed_trajectories = zeros(length(tv), 1);

if ~iscell(em)
    meanVals = mean(mean(data_array, 1), 3); % a 1 by # measured variables array
    wt = sum(meanVals)./meanVals; %hight mean = lower wt
    relWt = wt/sum(wt); % note that relWt is a row vector. 
    
    CONC_temp = zeros(length(tv), length(ms), size(dv,2));
    for i = 1:size(dv,2)
        sd = simulate(em, [exp(logp); dv(:,i)]);
        sd = resample(sd, tv);
        for j = 1:length(ms)
            measuredspecies = ms{j};
            spSD = selectbyname(sd, measuredspecies);
            summed_trajectories = sum(spSD.Data, 2);
            CONC_temp(:, j, i) = summed_trajectories;
        end
    end
    
    ExpData = data_array;
    if isequal(size(CONC_temp), size(ExpData)) % both must be #timepoints x
        % # measured species x # dose combinations
        relWt_tiled = repmat(relWt, size(CONC_temp,1), 1,  size(CONC_temp,3));
        residuals = relWt_tiled.*(CONC_temp - ExpData);
        res = residuals(:);
%         scaledsimdata = relWt_tiled.*(CONC_temp);
%         scaleddata = relWt_tiled.*(ExpData);
%         figure
%         for ccol = 1:2
%             for rrow = 1:4
%                 subplot(4, 2, (rrow-1)*2+ccol)
%                 plot(tv, CONC_temp(:,ccol, rrow), 'r', tv, scaledsimdata(:,ccol, rrow), 'b');
%                 hold on
%                 plot(tv, ExpData(:,ccol, rrow), ':r', tv, scaleddata(:,ccol, rrow), ':b');
%                 hold on 
%                 plot(tv, residuals(:,ccol, rrow), 'k');
%             end
%         end
        
%         figure
%         for ccol = 1:2
%             for rrow = 1:4
%                 subplot(4, 2, (rrow-1)*2+ccol)
%                 plot(tv, scaledsimdata(:,ccol, rrow), 'b');
%                 hold on
%                 plot(tv, scaleddata(:,ccol, rrow), ':b');
%                 hold on 
%                 plot(tv, residuals(:,ccol, rrow), 'k');
%             end
%         end
        
        % scale the residuals by relative weight. THis emphasizes
        % the species which are low conc,
        % and deemphasizes species which are high in magnitude
        % to make them all equally important
    else
        error('txtl_toolbox:objFun_multiDoseSpecies:incompatibleArrays',...
            'The simulation and experimental data need to be the same sizes.')
    end
else
    % start combined optimization mode % write this properly later
    % check that all the arrays are in the right format.
    % I think using cell arrays here is probably very slow. But lets try it
    % nonetheless for generality
    if ~iscell(data_array) || ~iscell(tv) || ~iscell(dv) ...
            || ~iscell(ms)
        error('txtl_toolbox:genresiduals:inputsIncompatible1',...
            'Cell array expected for data, time vector, dose values, and measured species')
    end
    assert(isequal(length(em), length(data_array), ...
        length(tv), length(dv), length(ms)), ...
        'txtl_toolbox:genresiduals:inputsIncompatible2',...
        'The number of cells in the input arrays must be equal')
    
    nOpt = length(em);
%     preallocate redidual array for speed

    nres = zeros(nOpt, 1);


    for kk = 1:nOpt
        mS = ms{kk};
        tv = tv{kk};
        nms = length(mS);
        nres(kk) = length(tv)*nms;
        
    end
    res = zeros(sum(nres), 1);
    
    for kk = 1:nOpt
        [~, ~, V] = find(multiparam(kk, :));
        lp = logp(V);
        da = data_array{kk};
        dIV = dv{kk};
        mS = ms{kk};
        tv = tv{kk};
        eMO = em{kk};
        
        meanVals = mean(da);
        % meanVals = mean(reshape(data_array, [size(data_array, 1) size(data_array, 3)]),1);
        wt = sum(meanVals)./meanVals; %hight mean = lower wt
        %!TODO is the mean the right statistic, or is the median or max a better statistic?
        relWt = wt/sum(wt);
        CONC_temp = zeros(length(tv)*size(dIV,1), length(mS));
        for i = 1:size(dIV,1)
            
            sd = simulate(eMO, [exp(lp); dIV(i,:)']);
            sd = resample(sd, tv);
            
            spSD = selectbyname(sd, mS);
            % !TODO :test if i can just feed in the cell of strings.
            CONC_temp((i-1)*length(tv) + 1:(i-1)*length(tv) +length(tv), :) ...
                = spSD.Data; % and then just do this.
        end
        ExpData = da;
        if isequal(size(CONC_temp), size(ExpData))
            residuals = repmat(relWt, size(CONC_temp,1),1).*(CONC_temp - ExpData);
            residx = (sum(nres(1:(kk-1)))+1):sum(nres(1:kk));
            res(residx) = residuals(:);
            % scale the residuals by relative weight. THis emphasizes
            % the species which are low conc,
            % and deemphasizes species which are high in magnitude
            % to make them all equally important
        else 
            error('simulation and experimental data must have same sizes')
        end
        
    end
    % can try interp1 or sd.resample to see which is faster.
    % I think deciding which strategy is faster will require some testing.
    
    % for now lets just use the naive sd.resample method.
    % If its abysmally slow, can try to modify this code to have interp1,
    % with the measured species indices predetermined.
end

