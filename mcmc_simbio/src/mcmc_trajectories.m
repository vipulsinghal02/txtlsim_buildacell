function [fighandle, varargout] = mcmc_trajectories(em, di, mi, marray,...
    titl, lgds, varargin)
%
% Plot the data time course trajectories and simulated model trajectories
% for each dose and measured species. The doses are arranged in rows of a
% subplot, or all collapsed into a single row. Each measured species has its
% own figure.
%
% Subplot arrangement:
% - Subplot Column: The column corresponds to a measured species. If there
% 		are multiple measured species, then multiple figures are generated.
% 		If the separateExpSim optional input parameter's value is set to
% 		true, then there are two columns, the first one corresponding the
% 		experimental data and the second corresponding to the simulated data.
% - Subplot Rows: Default is to use one row for each dose. But if the
% 		collateDoses input parameter is set to true, then all the doses
% 		get plotted on a single row.
%
% The experimental data and the simulation may be plotted in a few different
% ways: Mean + standard deviation, (curvewise) median + rest of the curves,
% just the curves, just the mean, just the median.
% The default is to use the median for the experimental data, and the mean
% + standard deviation for the simulated curves.
% The standard deviation is plotted as a shaded region. When the individual
% curves are plotted, they are plotted as thin lines: solid for experimental
% data and dotted for simulations.
%
% OUTPUTS: fighandle and tw optional outputs:
% varargout{1} = data array of dimensions nTimepoints x nMeasuredSpecies x
% nSimCurves x nDoseCombinations
% varargout{2} = idxnotused. This is the set of indices of the third
% (nSimCurves) dimension that have at least one dose that led to an error
% during the simulation. The corresponding vector in the da array will have
% all NaNs. 
%
% INPUTS
% em: The exported simbiology model object that is to be simulated.
% This is created when the mcmc_runsim is run, and the parameters, species an
% doses that can be set are predefined.
%
% data_info: This is the data info struct that contains the data and the
% related information.
%
% mcmc_info: This is the mcmc info struct that contains info on the mcmc
% run (including which species are dosed, measured etc.)
%
% marray: A numerical array of dimensions nParam x nWalkers x nSamp
% or nPoints x nParam.
%
% OPTIONAL NAME VALUE PAIRS
% collateDoses: Default is false. If true, all the doses get plotted on the
% same row of the subplots.
%
% separateExpSim: Default is false. If true, the experimental data and the
% simulation graphs are each given their own column.
%
% ExpMode: How to plot the experimental data. Can be 'mean', 'median'
% 'meanstd','medianstd', 'mediancurves', 'meancurves', 'curves' or 'none'.
% Default is 'median'.
%
% SimMode: How to plot the simulated data. Can be 'mean', 'median'
% 'meanstd','medianstd', 'mediancurves', 'meancurves', 'curves' or 'none'.
% Default is 'median'.
% If both the Exp and Sim Mode are none, or if the separateExpSim is set to
% true, and either mode are none, then an error is thrown.
%
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% How to specify TITLES and LEGENDS.
% Titles and Legends depend on the combination of CollatedDoses and separateExpSim
%
% CollatedDoses 1; separateExpSim 1
% The subplots within each species' figure are arranged as 1 x 2.
%
% Title: This should we a cell array of strings of size 1 x 2 x nMS, where nMS is the
% 		number of measured species. The first column corresponds to the experimental
% 		data, and the second to the simulated data. Each string should specify the
% 		measures species that is being displayed in the plot, and whether it is
% 		experimental or simulated data. Furthermore, You should also specify
% 		what kind of statistic is being shown: mean or median with standard
%		deviation or just sample curves etc.
% Legends: This should be a cell array of strings of dimension 1 x nICs or
% 		2 x nICs (the two rows corresponding to experimental data column,
% 		and simulated data column respectively). If either ExpMode or SimMode are
% 		'none', then the dimension must be 1 x nICs. If the dimension is 1 x nICs,
% 		and both ExpMode and SimMode are specified, then the same legend is used for
% 		both.
%
% CollatedDoses 1; separateExpSim 0
% Here the subplots are arranged as 1 by 1 for each species.
%
% Title: This should be a cell array of strings of size 1 x nMS, where nMS is the
% 		number of measured species. Each string corresponds to one measured species.
% 		You should also probably use this place to specify what the statistics being
% 		plotted are: (mean median or none) with (std, curves or none) for both
% 		experimental data and simulated data.
%
% Legends: This should be a cell array of strings of dimension 1 x (nICs + 1)
% 		The strings should describe the curves corresponding to Dose 1 of
% 		experimental data, Dose 1 of simulated data, Doses 2 to last of
% 		experimental data. For the simulated data doses 2 to last, the fact that
% 		the same colors are used for the same dose for both the simulated and
% 		experimental data should be used to interpret the simulated data curves.
% 		If either ExpMode or SimMode are 'none', then the dimension must
% 		be 1 x nICs. The same legend sting array is used for every measured
% 		species.
%
% CollatedDoses 0; separateExpSim 1
% Here the number of subplots is nICs x 2 per figure, and there is one
% figure for each species.
%
% Title: This should be a cell array of size nICs x 2 x nMS, where nMS
% 		is the number of measured species.
% Legends:  No legend array is needed.
%
% CollatedDoses 0; separateExpSim 0
% Number of subplots: nICs x 1 per figure, and there is one figure for
% each species.
%
% Title: This should be a cell array of size nICs x 1 x nMS, where nMS
% 		is the number of measured species.
% Legends: No legend array is needed. In the first dose, the experimental
% 		data and simulated data plots are labeled as 'exp' and 'sim'
% 		respectively.
%
%
%
%
%
%
% If SimMode of ExpMode are 'curves', then the curves are used in the legend.
% If Either of them are none, then no legend is used.
%
% TITLES:
%
%
%
% If either the ExpMode or the SimMode are 'curves', then the
%

% Specify the title and legend inputs as follows:
% For
% The legend that gets specified, if one is to be specified
% (see 'Conditions under which legends are included in the plots'),
% for both ExpMode and SimMode is as follows:
% 'mean', 'median': Specify legend for the mean or median curves.
%

%
% If collateDoses is true, then the legends are the doses, otherwise the
% titles have the dose information, and the legends depend on the summary
% and spread statistics displayed in each subplot.
%
% The plot titles are taken from the data info struct's measuredNames field.
% If that field is empty, then the measured species field of the mcmc info
% struct is used. Dosing information is taken from the dosedNames and the
% dosedVals fields of the data info struct, and if these are not populated,
% then it is taken from the dosedNames and dosedVals fields of the mcmc_info
% struct.
%
% When the collateDoses option is false, i.e., each dose is
% is plotted in a separate row, then the dose values are also used in the
% title string. If collateDoses is true, then the dose info is used to produce
% the legend strings. Depending on the values of the options ExpSummary, ExpSpread,
% SimSummary, and SimSpread, we also include legends for the corresponding lines
% (if there are any) for the first dose.
%
%
%
% --------------------------------------------------------------------------

% Copyright (c) 2018, Vipul Singhal, Caltech
% Permission is hereby granted, free of charge, to any person obtaining a copy
% of this software and associated documentation files (the "Software"), to deal
% in the Software without restriction, including without limitation the rights
% to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
% copies of the Software, and to permit persons to whom the Software is
% furnished to do so, subject to the following conditions:

% The above copyright notice and this permission notice shall be included in all
% copies or substantial portions of the Software.

% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
% IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
% FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
% AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
% LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
% OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
% SOFTWARE.

%%  precompute a bunch of things to set the defaults for this function
% Number of simulation curves to plot

p = inputParser;
p.addParameter('nSimCurves', 50, @isnumeric);
p.addParameter('collateDoses', false, @islogical);
p.addParameter('separateExpSim', false, @islogical);
p.addParameter('just_data_info', false, @islogical);
p.addParameter('subplot_arrangement', [], @isnumeric); % [nrows ncols]. 
% Must have nrows*ncols = ndoses. only gets used if collate doses is false
% Modes for the next two inputs: 'mean', 'median' 'meanstd', 'medianstd'
% 'mediancurves', 'curves'
p.addParameter('ExpMode', 'median', @ischar);
p.addParameter('SimMode', 'meanstd', @ischar);
p.addParameter('title', {}, @iscell);
p.addParameter('legends', {}, @iscell);
p.addParameter('savematlabfig', false, @islogical); 
% if this is true, then projdir and tstamp must be specified.
p.addParameter('savejpeg', false, @islogical); 
% if this is true, then projdir and tstamp must be specified.
p.addParameter('projdir', [], @ischar);
p.addParameter('tstamp', [], @ischar);
p.addParameter('extrafignamestring', [], @ischar);
p.parse(varargin{:})
p = p.Results;

% get the screen size for plotting. 
scrsz = get(0,'screensize');

if p.separateExpSim && (strcmp(p.ExpMode,'none') || strcmp(p.SimMode, 'none'))
    error(['Cant have unspecified Experimental data or Simulation Data '...
        'if the separateExpSim is set to true'])
end

if p.just_data_info
    % just plot the data in data info
    % One plot per measured species in each data info. All doses collated
    % in the same plot.
    fighandle = cell(length(di), 1);
    
    for dID = 1:length(di)
        currdi = di(dID);
        [expsummst, expspreadst] = computeDataStats(currdi.dataArray, p.ExpMode);
        %         dimensionLabels = currdi.dimensionLabels;
        %         expmax = computeMaxes(expsummst, ...
        % expspreadst, p.ExpMode, dimensionLabels);
        dNames = currdi.dosedNames;
        dVals = currdi.dosedVals;
        [ndNames, nICs] = size(dVals);
        assert(length(dNames) == ndNames);
        linehandle = zeros(nICs, 1);
        %         ptchhandle = zeros(nICs, 1);
        nMS = length(currdi.measuredNames);
        
        tv = currdi.timeVector;
        timeUnits = currdi.timeUnits;
        tv = converttosec(tv, timeUnits);
        colorz = parula(nICs+2);
        for msnum = 1:nMS
            fighandle{dID}(msnum) = figure( 'position',[50 50 scrsz(3)/1.7 scrsz(4)/1.4]);
            ax = gca;
            legendentry = cell(ndNames, 1);
            for i=1:nICs
                [ax, linehandle(i)] = ...
                    plotintoaxis(ax, p.ExpMode,...
                    tv, expsummst, expspreadst, ...
                    i, msnum, ...
                    'LineColor', colorz(i, :),...
                    'LineStyle', '-',...
                    'LineWidth', 2,...
                    'SpreadColor', colorz(i, :),...
                    'SpreadLineStyle', '--',...
                    'SpreadLineWidth', 0.5,...
                    'FaceAlpha', 0.25); % removed out arg: (doesnt work yet) , ptchhandle
                legendentry{i} = [];
                
                for dnID = 1:ndNames-1
                    legendentry{i} = ...
                    [legendentry{i} dNames{dnID} ' = ' num2str(dVals(dnID, i)) ', '];
                end
                
                legendentry{i} = ...
                [legendentry{i} dNames{ndNames} ' = ' num2str(dVals(ndNames,i ))];
            end
            legend(linehandle, legendentry);
            %title(currdi.measuredNames{msnum}{1:end}) % what was I thinking here?
            title(currdi.measuredNames{msnum})
        end
        
        
    end
    
    
else
    
    % Compute the y axis limits for each measured species. All the y axis for a
    % given measured species are set to the max y axis limits.
    
    %% set the number of curves to simulate
    % number of walkers is the second dimension of the 3D version of the
    % paraemter array OR if the parameter array is 2D, it is taken from the
    % mi array.
    if ndims(marray) == 3
        % compute the number of walkers
        nWalkers = size(marray, 2);
        m = marray(:,:)'; % the parameter array is now #points x #params
    elseif ismatrix(marray)
        % 		nWalkers = mi.nW;
        m = marray; % assume that the 2 dims are correct - npoints x nparams
    end
    
    % Number of curves to simulate is the minimum of the specified number and the
    % number of available walkers.
    % 	if p.nSimCurves > nWalkers
    % 		p.nSimCurves = nWalkers;
    %     end
    
    if p.nSimCurves > size(m, 1)
        p.nSimCurves = size(m, 1);
    end
    
    %% initialize things
    [nDSP, nICs] = size(mi.dosedVals); % number of dosed species,
    % and number of dose combinations
    
    dn = mi.dosedNames; %cell array
    dose  = mi.dosedVals';
    % recall from the data_info documentation:
    % 'dosedVals': A matrix of dose values of size
    % # of dosed species by # of dose combinations
    % thus dose has dimensions #combos x #species
    
    % convert time vector to seconds.
    tv = di.timeVector;
    timeUnits = di.timeUnits;
    tv = converttosec(tv, timeUnits);
    
    nts = length(tv);
    nMS = length(mi.measuredSpecies);
    
    % initialize arrays (nominal array is:
    % timepoints x measured outputs x nSimCurves x doses )
    da = zeros(nts, nMS, p.nSimCurves, nICs);
    ms = mi.measuredSpecies;
    dv = mi.dosedVals;
    %%
    % figure out how to incorporate Extract differences too.
    % !TODO
    
    %% Simulate the model for the given parameters.
    % set parameters in the model
    % for each dose simulate the model
    % simulate the model
    
    [da, idxnotused] = simulatecurves(em,m, p.nSimCurves, dose, tv, ms);
    % compute the simulation maxes, with the non - integrable points removed.
    pointstouse = setdiff(1:(p.nSimCurves),idxnotused);
    if isempty(pointstouse)
        error('none of the points are integrable. Something has gone wrong...')
    end
    %% Compute relevant statistics depending on what is specified in the inputs.
    % for the experimental data.
    [expsummst, expspreadst] = computeDataStats(di.dataArray, p.ExpMode);
    [simsummst, simspreadst] = computeDataStats(da(:,:,pointstouse, :), p.SimMode);
    
    dimensionLabels = {'time points', 'measured species', 'replicates', 'doses'};
    expmax = computeMaxes(expsummst, expspreadst, p.ExpMode, dimensionLabels);
    simmax = computeMaxes(simsummst, simspreadst, p.SimMode, dimensionLabels);
    expsimmax = max(expmax, simmax); %%
    colorz = parula(nICs+2);
    colorz2 = summer(nICs+2);
    colorz3 = winter(nICs+2);
    fighandle = zeros(nMS, 1);
    % 	titl = p.title;
    % 	legs = p.legends;
    
    
    
    
    for msnum = 1:nMS
        fighandle(msnum) = figure( 'position',[50 50 scrsz(3)/1.7 scrsz(4)/1.4]);;
        ax = gca;
        if p.collateDoses
            if p.separateExpSim
                % SEPARATEEXPSIM 1 , COLLATEDOSES 1
                ax = subplot(1, 2, 1); % get axis
                linehandle = zeros(nICs, 1);
                ptchhandle = zeros(nICs, 1);
                for i=1:nICs
                    [ax, linehandle(i)] = ...
                        plotintoaxis(ax, p.ExpMode,...
                        tv, expsummst, expspreadst, ...
                        i, msnum, ...
                        'LineColor', colorz(i, :),...
                        'LineStyle', '-',...
                        'LineWidth', 2,...
                        'SpreadColor', colorz(i, :),...
                        'SpreadLineStyle', '--',...
                        'SpreadLineWidth', 0.5,...
                        'FaceAlpha', 0.25);
                end
                
%                 legends([linehandle], legs(1,:));
                title(titl{1, 1, msnum});
                % !! p.legends must be a cell array of size 2 x nICs
                % p.title must be a 2 x nMS cell array containing titles
                % specifying experiment and simulation for a given measures species.
                % can talk about what the summary and spread statistic refer to.
                ax = subplot(1, 2, 2);
                for i=1:nICs
                    [ax, linehandle] = ...
                        plotintoaxis(ax, p.SimMode,...
                        tv, simsummst, simspreadst, ...
                        i, msnum, ...
                        'LineColor', colorz(i, :),...
                        'LineStyle', '--',...
                        'LineWidth', 2,...
                        'SpreadColor', colorz(i, :),...
                        'SpreadLineStyle', ':',...
                        'SpreadLineWidth', 0.5,...
                        'FaceAlpha', 0.25); % removed out arg: (doesnt work yet) , ptchhandle
                end
                % legends([linehandle(1); ptchhandle(1); linehandle(2:end)],
                % legs(2,:));
%                 legends([linehandle], legs(2,:));
                title(titl{1, 2, msnum})
                
            else
                % SEPARATEEXPSIM 0 , COLLATEDOSES 1
                
                linehandle = zeros(nICs, 1);
                ptchhandle = zeros(nICs, 1);
                
                linehandle2 = zeros(nICs, 1);
                ptchhandle2 = zeros(nICs, 1);
                % plot experimental data
                for i=1:nICs
                    [ax, linehandle] = ...
                        plotintoaxis(ax, p.ExpMode,...
                        tv, expsummst, expspreadst, ...
                        i, msnum, ...
                        'LineColor', colorz(i, :),...
                        'LineStyle', '-',...
                        'LineWidth', 2,...
                        'SpreadColor', colorz(i, :),...
                        'SpreadLineStyle', '--',...
                        'SpreadLineWidth', 0.5,...
                        'FaceAlpha', 0.25); % removed out arg: (doesnt work yet) , ptchhandle
                end
                % plot simulation data
                for i=1:nICs
                    [ax, linehandle2] = ...
                        plotintoaxis(ax, p.SimMode,...
                        tv, simsummst, simspreadst, ...
                        i, msnum, ...
                        'LineColor', colorz2(i, :),...
                        'LineStyle', '-.',...
                        'LineWidth', 1,...
                        'SpreadColor', colorz2(i, :),...
                        'SpreadLineStyle', ':',...
                        'SpreadLineWidth', 0.25,...
                        'FaceAlpha', 0.25); % removed out arg: (doesnt work yet) , ptchhandle2
                end
                title(titl{msnum});
                % !! p.legends must be a legend array saying exp dose 1,
                % sim dose 1, exp dose 2 to end.
                legends([linehandle(1); linehandle2(1); linehandle(2:end)],...
                    p.legends);
            end
        else
            if ~isempty(p.subplot_arrangement)
                assert(numel(p.subplot_arrangement) == 2);
                assert(prod(p.subplot_arrangement) == nICs);
                pCell = num2cell(p.subplot_arrangement);
                [nrows, ncols] = pCell{:};
                
                
                linehandle = zeros(nrows, ncols);
                ptchhandle = zeros(nrows, ncols);

                linehandle2 = zeros(nrows, ncols);
                ptchhandle2 = zeros(nrows, ncols);

                for i = 1:nICs
                    rowIX = floor((i-1)/ncols)+1;
                    colIX = i-floor((i-1)/ncols)*ncols;
                   
                    if p.separateExpSim
                        warning('Both subplot arrangement and separate experiment and sim are specified.\n ' ...
                            'Not going to do anything')
                        % future version : just make two figures. 
                        
                        % SEPARATEEXPSIM 1 , COLLATEDOSES 0,
                        % subplot_arrangement specified. 
%                         
%                         ind = (i-1)*2+1;
%                         ax = subplot(nICs, 2, ind); % exp data plot.
%                         [ax, linehandle] = ...
%                             plotintoaxis(ax, p.ExpMode,...
%                             tv, expsummst, expspreadst, ...
%                             i, msnum, ...
%                             'LineColor', colorz(i, :),...
%                             'LineStyle', '-',...
%                             'LineWidth', 2,...
%                             'SpreadColor', colorz(i, :),...
%                             'SpreadLineStyle', '--',...
%                             'SpreadLineWidth', 0.5,...
%                             'FaceAlpha', 0.25); % removed out arg: (doesnt work yet) , ptchhandle
%                         % set title
%                         title(titl{i, 1, msnum});
%                         ind = (i-1)*2+2;
%                         ax = subplot(nICs, 2, ind); % simulated data plot
%                         [ax, linehandle] = ...
%                             plotintoaxis(ax, p.SimMode,...
%                             tv, simsummst, simspreadst, ...
%                             i, msnum, ...
%                             'LineColor', colorz(i, :),...
%                             'LineStyle', '-.',...
%                             'LineWidth', 2,...
%                             'SpreadColor', colorz(i, :),...
%                             'SpreadLineStyle', ':',...
%                             'SpreadLineWidth', 0.5,...
%                             'FaceAlpha', 0.25); % removed out arg: (doesnt work yet) , ptchhandle
%                         % set title
%                         title(titl{i, 2, msnum});
%                         % in this case, the title has to be a cell array
%                         % of dimensions nICs x 2 x nMS
                    else
                        % % SEPARATEEXPSIM 0 , COLLATEDOSES 0
                        
                        
                        ax = subplot(nrows, ncols, i); % start plotting row first. ugh. why matlab, why?
                        [ax, linehandle(rowIX, colIX)] = plotintoaxis(ax, p.ExpMode,...
                            tv, expsummst, expspreadst, ...
                            i, msnum, ...
                            'LineColor', colorz(i, :),...
                            'LineStyle', '-',...
                            'LineWidth', 2,...
                            'SpreadColor', colorz(i, :),...
                            'SpreadLineStyle', '--',...
                            'SpreadLineWidth', 0.5,...
                            'FaceAlpha', 0.25); % removed out arg: (doesnt work yet) , , ptchhandle(i)
                        hold on
                        [ax, linehandle2(rowIX, colIX)] = ...
                            plotintoaxis(ax, p.SimMode,...
                            tv, simsummst, simspreadst, ...
                            i, msnum, ...
                            'LineColor', colorz2(i, :),...
                            'LineStyle', '-.',...
                            'LineWidth', 1,...
                            'SpreadColor', colorz2(i, :),...
                            'SpreadLineStyle', '--',...
                            'SpreadLineWidth', 0.25,...
                            'FaceAlpha', 0.25); % removed out arg: (doesnt work yet) ,, ptchhandle2(i)
                        hold on
                        % set titles
                        title(titl{i, 1, msnum});
                        % in this case, the p.title has to be a cell array
                        % of dimensions
                        % nICs x 1 x nMS
                        % set legends only for the top subplot.
                        % and distinguish the experimental and simulated data.
                        
                        if i == 1
                            legend([linehandle(1); linehandle2(1)],...
                                {'exp', 'sim'}, 'Location', 'SouthEast')
                        end
                    end
                end
                
            else
                for i = 1:nICs
                    if p.separateExpSim
                        
                        % SEPARATEEXPSIM 1 , COLLATEDOSES 0
                        
                        ind = (i-1)*2+1;
                        ax = subplot(nICs, 2, ind); % exp data plot.
                        [ax, linehandle] = ...
                            plotintoaxis(ax, p.ExpMode,...
                            tv/60, expsummst, expspreadst, ...
                            i, msnum, ...
                            'LineColor', colorz(i, :),...
                            'LineStyle', '-',...
                            'LineWidth', 2,...
                            'SpreadColor', colorz(i, :),...
                            'SpreadLineStyle', '--',...
                            'SpreadLineWidth', 0.5,...
                            'FaceAlpha', 0.25); % removed out arg: (doesnt work yet) , ptchhandle
                        ax.YLim = [0 expsimmax(msnum)];
                        % set title
                        title(titl{i, 1, msnum});
                        ind = (i-1)*2+2;
                        ax = subplot(nICs, 2, ind); % simulated data plot
                        [ax, linehandle] = ...
                            plotintoaxis(ax, p.SimMode,...
                            tv/60, simsummst, simspreadst, ...
                            i, msnum, ...
                            'LineColor', colorz(i, :),...
                            'LineStyle', '-.',...
                            'LineWidth', 2,...
                            'SpreadColor', colorz(i, :),...
                            'SpreadLineStyle', ':',...
                            'SpreadLineWidth', 0.5,...
                            'FaceAlpha', 0.25); % removed out arg: (doesnt work yet) , ptchhandle
                        ax.YLim = [0 expsimmax(msnum)];
                        % set title
                        title(titl{i, 2, msnum});
                        % in this case, the title has to be a cell array
                        % of dimensions nICs x 2 x nMS
                    else
                        % % SEPARATEEXPSIM 0 , COLLATEDOSES 0
                        
                        linehandle = zeros(nICs, 1);
                        ptchhandle = zeros(nICs, 1);
                        
                        linehandle2 = zeros(nICs, 1);
                        ptchhandle2 = zeros(nICs, 1);
                        
                        ax = subplot(nICs, 1, i);
                        [ax, linehandle(i)] = plotintoaxis(ax, p.ExpMode,...
                            tv/60, expsummst, expspreadst, ...
                            i, msnum, ...
                            'LineColor', colorz(i, :),...
                            'LineStyle', '-',...
                            'LineWidth', 2,...
                            'SpreadColor', colorz(i, :),...
                            'SpreadLineStyle', '--',...
                            'SpreadLineWidth', 0.5,...
                            'FaceAlpha', 0.25); % removed out arg: (doesnt work yet) , , ptchhandle(i)
                        hold on
                        [ax, linehandle2(i)] = ...
                            plotintoaxis(ax, p.SimMode,...
                            tv/60, simsummst, simspreadst, ...
                            i, msnum, ...
                            'LineColor', colorz2(i, :),...
                            'LineStyle', '-.',...
                            'LineWidth', 1,...
                            'SpreadColor', colorz2(i, :),...
                            'SpreadLineStyle', '--',...
                            'SpreadLineWidth', 0.25,...
                            'FaceAlpha', 0.25); % removed out arg: (doesnt work yet) ,, ptchhandle2(i)
                        ax.YLim = [0 expsimmax(msnum)];
                        hold on
                        % set titles
                        title(titl{i, 1, msnum});
                        % in this case, the p.title has to be a cell array
                        % of dimensions
                        % nICs x 1 x nMS
                        % set legends only for the top subplot.
                        % and distinguish the experimental and simulated data.
                        
                        if i == 1
                            legend([linehandle(1); linehandle2(1)],...
                                {'exp', 'sim'})
                        end
                    end
                end
                
            end
            
            
        end
        
        % save figure here
        if p.savematlabfig
            if isempty(p.projdir) || isempty(p.tstamp)
                warning('timestamp and project directory not specified. Nothing will be saved.')
            else
                specificproj = [p.projdir '/simdata_' p.tstamp];
                saveas(gcf, [specificproj '/traj' p.tstamp num2str(msnum) p.extrafignamestring]);
                
            end
        end
        
        if p.savejpeg
            if isempty(p.projdir) || isempty(p.tstamp)
                warning('timestamp and project directory not specified. Nothing will be saved.')
            else
                specificproj = [p.projdir '/simdata_' p.tstamp];
                print(gcf, '-djpeg', '-r200', [specificproj '/traj' p.tstamp num2str(msnum) p.extrafignamestring])
            end
        end
    end
end

if nargout == 2
    varargout{1} = da;
elseif nargout == 3
    varargout{1} = da;
    varargout{2} = idxnotused;
end





end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%                 %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%% LOCAL FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%                 %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function tv = converttosec(tv, timeUnits)
% convert weeks, days, hours, or minutes to seconds.
switch timeUnits
    case 'weeks'
        tv = tv*7*24*3600;
    case 'days'
        tv = tv*24*3600;
    case 'hours'
        tv = tv*3600;
    case 'minutes'
        tv = tv*60;
    case 'seconds'
        tv = tv;
end
end

function [timeDim, measuredDim, doseDim, replicateDim] =...
    dimensionLabelMaps(dimensionLabels)

% dimension lebels have to be the strings
% 'replicates', 'time points', 'measured species', 'doses'.
replicateDim = strcmp(dimensionLabels, 'replicates');
replicateDim = find(replicateDim);
if isempty(replicateDim)
    error('There is no ''replicates'' entry in the dimension labels array.')
end
timeDim = strcmp(dimensionLabels, 'time points');
timeDim = find(timeDim);
if isempty(timeDim)
    error('There is no ''time points'' entry in the dimension labels array.')
end
measuredDim = strcmp(dimensionLabels, 'measured species');
measuredDim = find(measuredDim);
if isempty(measuredDim)
    error('There is no ''measured species'' entry in the dimension labels array.')
end
doseDim = strcmp(dimensionLabels, 'doses');
doseDim = find(doseDim);
if isempty(doseDim)
    error('There is no ''doses'' entry in the dimension labels array.')
end

end

function datamax = computeMaxes(summst, spreadst, dispmode, dimensionLabels)
[tD, mD, dD, rD] = dimensionLabelMaps(dimensionLabels);

switch dispmode
    case 'mean'
        % summst must be a time x measured species x 1 x doses array.
        datamax = max(max(summst, [],1), [], 4); % max over time and doses.
    case 'median'
        % summst must be a time x measured species x 1 x doses array.
        datamax = max(max(summst, [],1), [], 4);
    case 'meanstd'
        % summst must be a time x measured species x 1 x doses array.
        % spreadst must be a time x measured species x 1 x doses array.
        summ_spread_st = summst + spreadst;
        datamax = max(max(summ_spread_st, [],1), [], 4);
    case 'meancurves'
        datamax = max(max(max(spreadst, [],1), [], 4), [], 3);
    case 'medianstd'
        % summst must be a time x measured species x 1 x doses array.
        % spreadst must be a time x measured species x 1 x doses array.
        summ_spread_st = summst + spreadst;
        datamax = max(max(summ_spread_st, [],1), [], 4);
    case 'mediancurves'
        datamax = max(max(max(cat(3,spreadst, summst), [],1), [], 4), [], 3);
    case 'curves'
        datamax = max(max(max(spreadst, [],1), [], 4), [], 3);
    otherwise
        error(['Invalid data display mode. Must be one of: ''mean'','...
            ' ''median'' , ''meanstd'',''medianstd'', ''mediancurves'', ''curves''.'])
end
end


%!!




%!!

function nonMedianCurves = allButMedianCurve(dataArray, ix)
% data array must be time x measuredspecies x replicates x doses

% Get the median curves. Can this be done via pure vector indexing?
% !UNIMPORTANT but could be fun to think about sometime.
nonMedianCurves = zeros(size(dataArray, 1), size(dataArray, 2), ...
    size(dataArray, 3)-1, size(dataArray, 4));
ixx = 1:size(dataArray, 3);
for i = 1:size(dataArray, 2)
    for j = 1:size(dataArray, 4)
        ixxx = setdiff(ixx , ix(1,i, 1,j)); % remove median curve index
        for k = 1:size(nonMedianCurves, 3)
            nonMedianCurves(:,i,k,j) = dataArray(:, i, ixxx(k), j);
        end
    end
end
end


function [title, legend] = titlelegend(p, di, mi)
end

function [ax, varargout] = plotintoaxis(ax, mode, tv,...
    summst, spreadst, dosei, msi, varargin)
% plot a single summary statistic and spread statistic on a given axis
% optional name value pair arguments have names:
%
% 'LineColor'
% 'LineStyle'
% 'LineWidth'
% 'SpreadColor'
% 'SpreadLineStyle'
% 'SpreadLineWidth'
% 'FaceAlpha'
%
p = inputParser;
p.addParameter('LineColor', parula(1), @isnumeric);
p.addParameter('LineStyle','-' , @ischar);
p.addParameter('LineWidth', 2, @isnumeric);
p.addParameter('SpreadColor', summer(1), @isnumeric);
p.addParameter('SpreadLineStyle', ':', @ischar); %  gets used if the spread is curves.
p.addParameter('SpreadLineWidth', 2, @isnumeric);

% gets used if the spread is standard deviation
p.addParameter('FaceAlpha', 0.25, @isnumeric);

% shading.
p.parse(varargin{:})
p = p.Results;

% make the axes to be plotted on current, and bring the relevant figure into focus
axes(ax);

if strcmp(mode, 'meanstd') || strcmp(mode, 'medianstd')
    % plot summary and std as the spread.
    [linehandle, ptchhandle] =...
        boundedline(tv, summst(:,msi,1, dosei),...
        spreadst(:,msi,1, dosei));
    set(ptchhandle, ...
        'FaceColor', p.SpreadColor,...
        'FaceAlpha', p.FaceAlpha);
    set(linehandle, ...
        'Color', p.LineColor,...
        'LineStyle', p.LineStyle,...
        'LineWidth', p.LineWidth);
    hold on
    outstuff = {linehandle, ptchhandle};
    
elseif strcmp(mode, 'mean') || strcmp(mode, 'median')
    % plot only the mean or median (no spread)
    [linehandle] = plot(tv,	summst(:,msi,1, dosei));
    set(linehandle, ...
        'Color', p.LineColor,...
        'LineStyle', p.LineStyle,...
        'LineWidth', p.LineWidth);
    hold on
    outstuff = {linehandle};
    
elseif strcmp(mode, 'mediancurves') || strcmp(mode, 'meancurves')
    
    [linehandle] = plot(tv, ...
        summst(:,msi,1, dosei));
    set(linehandle, ...
        'Color', p.LineColor,...
        'LineStyle', p.LineStyle,...
        'LineWidth', p.LineWidth);
    hold on
    
    spreadhandles = zeros(size(spreadst, 3), 1);
    for j = 1:size(spreadst, 3) % number of replicates
        spreadhandles(j) = plot(tv, spreadst(:,msi,j, dosei));
        set(spreadhandles(j),...
            'Color', p.SpreadColor,...
            'LineStyle', p.SpreadLineStyle,...
            'LineWidth', p.SpreadLineWidth);
        hold on
    end
    outstuff = {linehandle, spreadhandles};
    
elseif strcmp(mode, 'curves')
    spreadhandles = zeros(size(spreadst, 3),1);
    for j = 1:size(spreadst, 3)
        spreadhandles(j) = plot(tv, ...
            spreadst(:,msi,j, dosei));
        set(spreadhandles(j),...
            'Color', p.SpreadColor,...
            'LineStyle', p.SpreadLineStyle,...
            'LineWidth', p.SpreadLineWidth);
        hold on
    end
    outstuff = {spreadhandles(1)};
else
    warning(['No valid exp data plotting format'...
        ' specified. No experimental data will be plotted.'])
end
hold on

% process variable outputs
nout = nargout-1;
varargout(1:nout) = outstuff(1:nout);



end













