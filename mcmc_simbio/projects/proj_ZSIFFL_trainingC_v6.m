function [mi,mai, ri, tstamp, projdir, di]  = proj_ZSIFFL_trainingC_v6(varargin)
% data collected by Zach Sun and Shaobin Guo. 
% Vipul Singhal,
% California Institute of Technology
% 2019

p = inputParser;
p.addParameter('prevtstamp', []);
p.addParameter('prevtstampID', []);
p.addParameter('stepSize', 1.4);
p.addParameter('nW', 30);
p.addParameter('nPoints', 30*100);
p.addParameter('thinning', 1);
p.addParameter('nIter', 2);
p.addParameter('pausemode', true);
p.addParameter('parallel', false);
p.addParameter('stdev', 1);
p.addParameter('poolsize', []);
p.addParameter('multiplier', 1);
p.addParameter('stepLadder', linspace(1.1, 1, 4), @isnumeric); % A vector of multipliers for the
% step size. Must have length > 0.5*nIter, since only the first nIter/2
% iterations get their step sizes changed.
% if stepLadder is specified, the multiplier is automatically set to 1.
p.addParameter('literalStepLadder', false)
p.addParameter('temperatureLadder', [0.00005]); % can be a boolean: true or false,
% or can be a vector of multipliers to allow for a simulated annealing type
% approach
p.parse(varargin{:});
p = p.Results;

% data_init
% proj_acs_dsg2014_regen_A('nW', 50, 'nPoints', 50*10*5, 'nIter', 5, 'parallel', false, 'multiplier', 2, 'thinning', 10)
% proj_acs_dsg2014_regen_A('nW', 6400, 'nPoints', 6400*10*20, 'nIter', 20, 'poolsize', 36, 'multiplier', 3, 'thinning', 10)
%% construct simbiology model object(s)
mlac = model_txtl_pLacdeGFP;
mlas = model_txtl_pLacLasR_pLasdeGFP;
mtet = model_txtl_ptetdeGFP_pLactetR_aTc;

%% setup the mcmc_info struct
mcmc_info = mcmc_info_ZSIFFL_training_fullC_v6(mtet, mlac, mlas);

mi = mcmc_info.model_info;

%% setup the data_info struct
di = data_ZSIFFL;
% modify di to only contain the mRNA data.
% di.dataArray = di.dataArray(:, 1, :, :); % pick out only the mrna
% di.measuredNames = di.measuredNames(1);
% di.dataUnits = di.dataUnits(1);
% di.dataInfo = ['Modified to only have mRNA data. \n',...
% di.dataInfo];

%     Run the MCMC
if ~isempty(p.stepSize)
    mcmc_info.runsim_info.stepSize = p.stepSize;
end
if ~isempty(p.nW)
    mcmc_info.runsim_info.nW = p.nW;
end
if ~isempty(p.nPoints)
    mcmc_info.runsim_info.nPoints = p.nPoints;
end
if ~isempty(p.thinning)
    mcmc_info.runsim_info.thinning = p.thinning;
end
if ~isempty(p.nIter)
    mcmc_info.runsim_info.nIter = p.nIter;
end
if ~isempty(p.parallel)
    mcmc_info.runsim_info.parallel = p.parallel;
end


if mcmc_info.runsim_info.parallel
    if ~isempty(p.poolsize)
        delete(gcp('nocreate'))
        parpool(p.poolsize)
    else
        delete(gcp('nocreate'))
        parpool
    end
end
if ~isempty(p.stdev)
    mcmc_info.runsim_info.stdev = p.stdev;
end
if ~isempty(p.multiplier)
    if ~isempty(p.stepLadder) && isvector(p.stepLadder) && isnumeric(p.stepLadder)
        p.multiplier = 1;
    end
    
end

ri = mcmc_info.runsim_info;
mai = mcmc_info.master_info;

%% Compute total signal strength
% compute total signal for the data of interest, taking dose and
% topology weighting into account.
tsig = 0; % total signal.
for kk = 1:length(mi)
    currda = di(mi(kk).dataToMapTo).dataArray; % current data array
    % renormalize the data to make the weight of the measured
    % species the same.
    meanVals = mean(mean(mean(currda, 1), 3), 4);
    wt = sum(meanVals)./meanVals; %hight mean = lower wt
    relWt = wt/sum(wt);
    for jj = 1:size(currda, 2)
        currda(:, jj, :, :) = relWt(jj)*currda(:, jj, :, :);
    end
    dv = mi(kk).dosedVals;
    % reweight data by dose and experient.
    if isfield(mi(kk), 'experimentWeighting')
        if ~isempty(mi(kk).experimentWeighting)
            if isfield(mi(kk), 'doseWeighting')
                if ~isempty(mi(kk).doseWeighting)
                    if isequal(size(mi(kk).doseWeighting), [1 size(dv,2)])
                        for ii = 1:size(dv,2)
                            currda(:, :, :, ii) = ...
                                mi(kk).doseWeighting(ii)...
                                *currda(:, :, :, ii)...
                                *mi(kk).experimentWeighting;
                        end
                    else
                        error('Doses weights must be the same length as the number of doses')
                    end
                end
            else
                currda = ...
                    currda(:, :, :, ii)...
                    *mi(kk).experimentWeighting;
            end
        end
    elseif isfield(mi(kk), 'doseWeighting')
        if ~isempty(mi(kk).doseWeighting)
            if isequal(size(mi(kk).doseWeighting), [1 size(dv,2)])
                for ii = 1:size(dv,2)
                    currda(:, :, :, ii) = ...
                        mi(kk).doseWeighting(ii)...
                        *currda(:, :, :, ii);
                end
            else
                error('Doses weights must be the same length as the number of doses')
            end
        end
    end
    tsig = tsig + sum(sum(sum(sum(currda)))); % total signal is "tsig"
end

%% run the mcmc simulations
% prevtstamp = {'20180120_172922'};
% simID = {'1'};
% marray = mcmc_get_walkers(prevtstamp, {simID}, projdir);
% mtemp = marray(:,:);







if islogical(p.temperatureLadder) % if p.temperatureLadder is logical
    if ~p.temperatureLadder % if p.temperatureLadder is false
        %% initialize the directory where things are stored.
        [tstamp, projdir, st] = project_init;
        if isempty(p.prevtstamp)
            mi = mcmc_runsim_v2(tstamp, projdir, di, mcmc_info,...
                'InitialDistribution', 'LHS', 'multiplier', p.multiplier,...
                'stepLadder', p.stepLadder,...
                'pausemode', p.pausemode,...
                'literalStepLadder', p.literalStepLadder);
        else
            specificprojdir = [projdir '\simdata_' p.prevtstamp];
            SS = load([specificprojdir '\full_variable_set_' p.prevtstamp], 'mcmc_info');
            if isempty(p.prevtstampID)
                marray = mcmc_get_walkers({p.prevtstamp}, {SS.mcmc_info.runsim_info.nIter},...
                    projdir);
            else
                marray = mcmc_get_walkers({p.prevtstamp}, {p.prevtstampID},...
                    projdir);
            end
            
            % assume the projdir where this data is stored is the same one as the
            % one created at the start of this file
            
            pID = 1:length(mai.estNames);
            marray_cut = mcmc_cut(marray, pID, flipud((mai.paramRanges)'));
            if size(marray_cut, 2) < ri.nW
                warning('too few initial points, using a few timesteps from previous runs to create initial walker positions.');
                walker_timepoints = ceil(linspace(ceil(size(marray_cut,3))/2, size(marray_cut,3), ceil(ri.nW/size(marray_cut, 2))))
                minit = marray_cut(:,:, walker_timepoints(1));
                for i = 2:length(walker_timepoints)
                    minit = [minit marray_cut(:,:,walker_timepoints(i)) ];
                end
                minit = minit(:, 1:ri.nW);
            else % there are enough points, just pick the number needed.
                minit = marray_cut(:,1:ri.nW,end);
            end
            % now run the simulation.
            
            mi = mcmc_runsim_v2(tstamp, projdir, di, mcmc_info,...
                'UserInitialize', minit, 'multiplier', 2,...
                'pausemode', p.pausemode, 'stepLadder', p.stepLadder,...
                'prevtstamp', p.prevtstamp,...
                'literalStepLadder', p.literalStepLadder);
        end
    else % if p.temperatureLadder is true
        % begin temperature ladder: Initial temp is 10% of the total signal
        % then we do 0.1%, then finally 0.001%
        disp('Using temperature ladder for MCMC at the following temperatures.');
        tladder = tsig*[0.001, 0.0001, 0.00001]
        percentLadder = {'_pt1pct', '_pt01pct', '_pt001pct'};
        tstamp = datestr(now, 'yyyymmdd_HHMMSS');
        for ll = 1:length(tladder)
            %         tLaddStr = num2str(tladder(ll));
            %         if ~isempty(regexp(tLaddStr, '\.', 'ONCE'))
            %             dotLoc = regexp(tLaddStr, '\.');
            %             tLaddStr(dotLoc) = 'p';
            %         end
            
            tstamp_appended = [tstamp '_' num2str(ll) '_' percentLadder{ll}];
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%% Set the standard deviation %%%
            mcmc_info.runsim_info.stdev = tladder(ll);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            [~, projdir, ~] = project_init('saveStr', tstamp_appended);
            if ll ==1 % first temperature. 
                if isempty(p.prevtstamp)
                    mi = mcmc_runsim_v2(tstamp_appended, projdir, di, mcmc_info,...
                        'InitialDistribution', 'LHS', 'multiplier', p.multiplier,...
                        'pausemode', p.pausemode,...
                        'stepLadder', p.stepLadder,...
                        'literalStepLadder', p.literalStepLadder);
                else
                    specificprojdir = [projdir '\simdata_' p.prevtstamp];
                    SS = load([specificprojdir '\full_variable_set_' p.prevtstamp], 'mcmc_info');
                    if isempty(p.prevtstampID)
                        marray = mcmc_get_walkers({p.prevtstamp}, {SS.mcmc_info.runsim_info.nIter},...
                            projdir);
                    else
                        marray = mcmc_get_walkers({p.prevtstamp}, {p.prevtstampID},...
                            projdir);
                    end
                    % assume the projdir where this data is stored is the same one as the
                    % one created at the start of this file
                    
                    pID = 1:length(mai.estNames);
                    marray_cut = mcmc_cut(marray, pID, flipud((mai.paramRanges)'));
                    if size(marray_cut, 2) < ri.nW
                        warning('too few initial points, using a few timesteps from previous runs to create initial walker positions.');
                        walker_timepoints = ceil(linspace(ceil(size(marray_cut,3))/2, size(marray_cut,3), ceil(ri.nW/size(marray_cut, 2))))
                        minit = marray_cut(:,:, walker_timepoints(1));
                        for i = 2:length(walker_timepoints)
                            minit = [minit marray_cut(:,:,walker_timepoints(i)) ];
                        end
                        minit = minit(:, 1:ri.nW);
                    else % there are enough points, just pick the number needed.
                        minit = marray_cut(:,1:ri.nW,end);
                    end
                    
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %%% refresh parallel pool to free memory on windows... (sigh)###
                    if mcmc_info.runsim_info.parallel
                        if ~isempty(p.poolsize)
                            delete(gcp('nocreate'))
                            parpool(p.poolsize)
                        else
                            delete(gcp('nocreate'))
                            parpool
                        end
                    end
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    
                    
                    % now run the simulation.
                    mi = mcmc_runsim_v2(tstamp_appended, projdir, di, mcmc_info,...
                        'UserInitialize', minit, 'multiplier', 1, ...
                        'pausemode', p.pausemode,...
                        'stepLadder', p.stepLadder, 'prevtstamp', p.prevtstamp,...
                        'literalStepLadder', p.literalStepLadder);
                end
            else % subsequent temperatures. 
                %                 prevtstamp = [percentLadder{ll-1} tstamp];
                specificprojdir = [projdir '\simdata_' prevtstamp];
                SS = load([specificprojdir '\full_variable_set_' prevtstamp], 'mcmc_info');
                
                marray = mcmc_get_walkers({p.prevtstamp}, {SS.mcmc_info.runsim_info.nIter},...
                    projdir);
                
                % assume the projdir where this data is stored is the same one as the
                % one created at the start of this file
                
                pID = 1:length(mai.estNames);
                marray_cut = mcmc_cut(marray, pID, flipud((mai.paramRanges)'));
                minit = marray_cut(:,1:ri.nW,end);
                
                
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%% refresh parallel pool to free memory on windows... (sigh)###
                if mcmc_info.runsim_info.parallel
                    if ~isempty(p.poolsize)
                        delete(gcp('nocreate'))
                        parpool(p.poolsize)
                    else
                        delete(gcp('nocreate'))
                        parpool
                    end
                end
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % now run the simulation.
                mi = mcmc_runsim_v2(tstamp_appended, projdir, di, mcmc_info,...
                    'UserInitialize', minit, 'multiplier', 1,...
                    'pausemode', p.pausemode,...
                    'stepLadder', p.stepLadder, 'prevtstamp', prevtstamp,...
                    'literalStepLadder', p.literalStepLadder);
            end
            % define tstamp for the next iter.
            prevtstamp = tstamp_appended;
            
        end
    end
    
    %p.temperatureLadder NOT logical
elseif isnumeric(p.temperatureLadder) && isvector(p.temperatureLadder)
    
    % begin temperature ladder:
    disp('Using temperature ladder for MCMC at the following temperatures.');
    tladder = tsig*p.temperatureLadder
    tstamp = datestr(now, 'yyyymmdd_HHMMSS');
    for ll = 1:length(tladder)
        if tladder(ll)>=1
            tLaddStr = num2str(round(tladder(ll)));
            fprintf('Current temperature: %d \n', round(tladder(ll)));
        else
            tLaddStr = num2str(tladder(ll));
            if ~isempty(regexp(tLaddStr, '\.', 'ONCE'))
                dotLoc = regexp(tLaddStr, '\.');
                tLaddStr(dotLoc) = 'p';
            end
            fprintf('Current temperature: %d \n', tladder(ll));
        end
        tstamp_appended = [tstamp '_' num2str(ll) '_' tLaddStr];
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% Set the standard deviation %%%
        mcmc_info.runsim_info.stdev = tladder(ll);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        [~, projdir, ~] = project_init('saveStr', tstamp_appended);
        if ll ==1
            % a previous time stamp is NOT provided, start with latin
            % hypercube sampling.
            if isempty(p.prevtstamp)
                mi = mcmc_runsim_v2(tstamp_appended, projdir, di, mcmc_info,...
                    'InitialDistribution', 'LHS', 'multiplier', p.multiplier,...
                    'stepLadder', p.stepLadder,...
                    'pausemode', p.pausemode,...
                    'literalStepLadder', p.literalStepLadder);
            else
                % a previous time stamp IS provided, use it as a starting
                % point. The data for that time stamp must be in the same
                % directory as the projdir for this project.
                specificprojdir = [projdir '\simdata_' p.prevtstamp];
                SS = load([specificprojdir '\full_variable_set_' p.prevtstamp], 'mcmc_info');
                if isempty(p.prevtstampID)
                    marray = mcmc_get_walkers({p.prevtstamp}, {SS.mcmc_info.runsim_info.nIter},...
                        projdir);
                else
                    marray = mcmc_get_walkers({p.prevtstamp}, {p.prevtstampID},...
                        projdir);
                end
                
                
%                 marray = mcmc_get_walkers({p.prevtstamp}, {SS.mcmc_info.runsim_info.nIter},...
%                     projdir);
                pID = 1:length(mai.estNames);
                marray_cut = mcmc_cut(marray, pID, flipud((mai.paramRanges)'));
                if size(marray_cut, 2) < ri.nW
                    warning('too few initial points, using a few timesteps from previous runs to create initial walker positions.');
                    walker_timepoints = ceil(linspace(ceil(size(marray_cut,3))/2, size(marray_cut,3), ceil(ri.nW/size(marray_cut, 2))))
                    minit = marray_cut(:,:, walker_timepoints(1));
                    for i = 2:length(walker_timepoints)
                        minit = [minit marray_cut(:,:,walker_timepoints(i)) ];
                    end
                    minit = minit(:, 1:ri.nW);
                else % there are enough points, just pick the number needed.
                    minit = marray_cut(:,1:ri.nW,end);
                end
                % now run the simulation.
                mi = mcmc_runsim_v2(tstamp_appended, projdir, di, mcmc_info,...
                    'UserInitialize', minit, 'multiplier', 1,...
                    'pausemode', p.pausemode,...
                    'stepLadder', p.stepLadder, 'prevtstamp', p.prevtstamp,...
                    'literalStepLadder', p.literalStepLadder);
            end
        else
            % get the final location of the walkers from the previous
            % iteration (ll - 1).
            specificprojdir = [projdir '\simdata_' prevtstamp];
            SS = load([specificprojdir '\full_variable_set_' prevtstamp], 'mcmc_info');
            marray = mcmc_get_walkers({prevtstamp}, {SS.mcmc_info.runsim_info.nIter},...
                projdir);
            % assume the projdir where this data is stored is the same one as the
            % one created at the start of this file

                    
            pID = 1:length(mai.estNames);
            marray_cut = mcmc_cut(marray, pID, flipud((mai.paramRanges)'));
            minit = marray_cut(:,1:ri.nW,end);
            % now run the simulation.
            mi = mcmc_runsim_v2(tstamp_appended, projdir, di, mcmc_info,...
                'UserInitialize', minit, 'multiplier', 1,...
                'pausemode', p.pausemode,...
                'stepLadder', p.stepLadder, 'prevtstamp', prevtstamp,...
                'literalStepLadder', p.literalStepLadder);
        end
        % define tstamp for the next iter.
        prevtstamp = tstamp_appended;
        
    end
    
    
else
    error('temterature ladder need to be specified as a boolean or a numeric vector')
end



% copy paste into the terminal if using matlab with no display.
% make sure these stay commented!!!
% proj_acs_dsg2014_regen_A('stepSize', 1.3, 'nW', 1000, 'nPoints', 100*1000,...
%     'thinning', 2, 'nIter', 2, 'parallel', true, 'poolsize', 4, 'temperatureLadder', true)
%
% proj_acs_dsg2014_regen_A('stepSize', 1.4, 'nW', 16*3200, 'nPoints', ...
%     16*3200*10*10, 'thinning', 1, 'nIter', 5, 'parallel', true, ...
%     'poolsize', 48, 'temperatureLadder', true)

% proj_acs_dsg2014_regen_A('stepSize', 1.2, 'nW', 800, 'nPoints', ...
%     800*10*20, 'thinning', 1, 'nIter', 10, 'parallel', true, ...
%     'poolsize', 23, 'temperatureLadder', [0.5, 0.005 0.05 0.0005 0.005 0.00005],...
%     'stepLadder', linspace(2, 1, 5))
%
% proj_acs_dsg2014_regen_A('stepSize', 1.2,...
%     'nW', 50,...
%     'nPoints', 10*3,...
%     'thinning', 1,...
%     'nIter', 4,...
%     'parallel', true, ...
%     'poolsize', 4,...
%     'temperatureLadder', [0.5, 0.005 0.05 0.0005],...
%     'stepLadder', linspace(2, 1, 2))

end

