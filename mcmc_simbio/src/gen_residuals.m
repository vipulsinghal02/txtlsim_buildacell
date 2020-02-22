function res = gen_residuals(logp, exportedMdlObj, data_array, timevec, ...
    dosedInitVals, measuredSpecies, varargin)
% generate residuals for MCMC
% dosedInitVals needs to be a nDoseCombinations x nSpeciesToDose matrix.
% the concentration of this matrix needs some work from the dosedArray struct.
% logp should be a column vector. in exportedMdlObj, we already have the things
% that can be varied to be just the parameters to be estimated and the species that can be dosed, in that order.
% measuresSpecies is a cell array of strings.
% (unlike the struct in objFun_multiDoseSpecies.m)
% data_array is a length(timevec) * nDoseCombinations by nMeasuredSPecies array.
% the order of the data is the same as the order of the dosing, and needs to be maintained carefully
%
p = inputParser;
% estimation structure is a matrix that specifies which parameters in logp
% to apply to the individual exported model objects. Ie, which parameters
% in logp get estimated for each (model, data, timevec, dose, measured
% species) tuple.
p.addParameter('multiopt_params', [], @isnumeric);
p.parse(varargin{:})
p = p.Results;
multiparam = p.multiopt_params; 

if ~iscell(exportedMdlObj)
    
    meanVals = mean(data_array);
    
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
        CONC_temp((i-1)*length(timevec) + 1:(i-1)*length(timevec) +...
            length(timevec), :) = spSD.Data; % and then just do this.
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
else
    % start combined optimization mode
    % check that all the arrays are in the right format.
    % I think using cell arrays here is probably very slow. But lets try it
    % nonetheless for generality
    if ~iscell(data_array) || ~iscell(timevec) || ~iscell(dosedInitVals) ...
            || ~iscell(measuredSpecies)
        error('txtl_toolbox:genresiduals:inputsIncompatible1',...
            'Cell array expected for data, time vector, dose values, and measured species')
    end
    assert(isequal(length(exportedMdlObj), length(data_array), ...
        length(timevec), length(dosedInitVals), length(measuredSpecies)), ...
        'txtl_toolbox:genresiduals:inputsIncompatible2',...
        'The number of cells in the input arrays must be equal')
    
    nOpt = length(exportedMdlObj);
%     preallocate redidual array for speed

nres = zeros(nOpt, 1);


    for kk = 1:nOpt
        mS = measuredSpecies{kk};
        tv = timevec{kk};
        nms = length(mS);
        nres(kk) = length(tv)*nms;
        
    end
    res = zeros(sum(nres), 1);
    
    for kk = 1:nOpt
        [~, ~, V] = find(multiparam(kk, :));
        lp = logp(V);
        da = data_array{kk};
        dIV = dosedInitVals{kk};
        mS = measuredSpecies{kk};
        tv = timevec{kk};
        eMO = exportedMdlObj{kk};
        
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
end

