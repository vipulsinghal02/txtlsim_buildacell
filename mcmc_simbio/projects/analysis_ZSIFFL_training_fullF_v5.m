

clear all
close all
saveFinalFigs = 'D:\Dropbox\Documents\a_Journal_Papers\Drafts\txtl_bmc_bioinformatics\figs\fullFv5\'

finafigmode = true
% Set working directory to the txtlsim toolbox directory.
projdir = [pwd '\mcmc_simbio\projects\proj_ZSIFFL_trainingF_v5'];
addpath(projdir)
sls = regexp(projdir, '\', 'split');
extrastring = sls{end};
jpgsave = true;
figsave = false;

% Load model, mcmc_info, and data_info.
% construct simbiology model object(s)
mtet = model_txtl_ptetdeGFP_pLactetR_aTc;
mlas = model_txtl_pLacLasR_pLasdeGFP;
mlac = model_txtl_pLacdeGFP;
% setup the mcmc_info struct
mcmc_info = mcmc_info_ZSIFFL_training_fullF_v5(mtet, mlac, mlas);
di = data_ZSIFFL;
mi = mcmc_info.model_info;
ri = mcmc_info.runsim_info;
mai = mcmc_info.master_info;
% plot data from existing simulations.
tsIDtouse = 1;
plotflag = true;
switch tsIDtouse
    case 1 % this is after including the termination parameters.
        ts1 = '20200528_165131_1_738';
        ts2 = '20200528_165131_2_554';
        ts3 = '20200528_165131_3_369';
        ts4 = '20200529_040611_1_369';
        ts5 = '20200529_040611_2_258';
        ts6 = '20200529_040611_3_148';        
        ts7 = '20200529_171720_1_148';
        ts8 = '20200529_171720_2_74';
        
        tstamp = {ts1 ts2 ts3 ts4 ts5 ts6 ts7 ts8};% 
        nIterID = {1:10 1:10 1 1:10 1:10 1:5 1:10 1:10};% 
        load([projdir '\simdata_' tstamp{end} '\full_variable_set_' tstamp{end} '.mat'], ...
            'mi',...
            'mcmc_info', 'data_info',  'ri');
end
tsToSave = tstamp{end};
mai.masterVector
marray_full = mcmc_get_walkers(tstamp,nIterID, projdir);
marray = marray_full(:,:,1:end);
clear marray_full
parnames = [...
    {'pol_{Kd, tet}'}
    {'pol_{Kd, lac}'}
    {'pol_{term}'}
    {'Ribo_{term}'}
    {'pol'}
    {'Ribo'}
    {'pol_{Kd,las}'}
    {'plas_{tf, Kd}'}
    {'plas-pol_{tf, Kd}'}    ];

mcmc_plot(marray(:, 1:end,1:end), parnames(:),...
    'savematlabfig', false, 'savejpeg', false,...
    'projdir', projdir);

%  [activeNames(estParamsIX,1) mat2cell(log(cell2mat(activeNames(estParamsIX,3))), ones(14,1), ones(1, 2))]
%
% ans =
%
%   14×3 cell array
%
%     {'TX_elong_glob'                  }    {[  0]}    {[ 5]}
%     {'TL_elong_glob'                  }    {[  0]}    {[ 6]}
%     {'TXTL_PTET_RNAPbound_Kd'         }    {[  0]}    {[22]}
%     {'TXTL_PLAC_RNAPbound_Kd'         }    {[  5]}    {[17]}
%     {'TXTL_RNAPBOUND_TERMINATION_RATE'}    {[  0]}    {[12]}
%     {'TXTL_RIBOBOUND_TERMINATION_RATE'}    {[  0]}    {[12]}
%     {'RNAP'                           }    {[  1]}    {[15]}
%     {'Ribo'                           }    {[  1]}    {[15]}
%     {'TXTL_PLAS_RNAPbound_Kd'         }    {[ 25]}    {[40]}
%     {'TXTL_PLAS_TFBIND_Kd'            }    {[  0]}    {[10]}
%     {'TXTL_PLAS_TFRNAPbound_Kd'       }    {[  0]}    {[15]}

% if plotflag
% % % %
% % % % close all
% % % % % %    close all
% % % % %%     % Plot trace and corner (posterior distribution) plots
%     mcmc_plot(marray(:, 1:end,1:4:end), parnames(:),...
%         'savematlabfig', figsave, 'savejpeg', jpgsave,...
%         'projdir', projdir, 'tstamp', tsToSave,...
%         'extrafignamestring', 'Without_transient');
% % % % %%
%
%%
mcmc_plot(marray(:, 1:end,end-300:50:end), parnames(:),...
    'plotChains', false, ...
    'savematlabfig', figsave, 'savejpeg', jpgsave,...
    'projdir', projdir, 'tstamp', tsToSave, 'extrafignamestring', 'Burned_in');
%%

close all

mcmc_plot(marray(:, 1:end,1:end), parnames(:),...
    'plotDistribution', false, ...
    'savematlabfig', figsave, 'savejpeg', jpgsave,...
    'projdir', projdir, 'tstamp', tsToSave, 'extrafignamestring', 'notsofewChains_all_iters');

% % % % %
%%
%
%
% parIDs = [4 5 11];
% %     parRanges(parIDs, :) = [-1 0 ;
% %         -4 0;
% %         11 15]
%     parRanges(parIDs, :) = [-1 0 ;
%         -3 -1;
%         11.5 14.5]
%
%     marray_cut = mcmc_cut(marray(:, 1:end,(end - 40):end), parIDs, flipud((parRanges(parIDs, :))'));
%
%     mcut0 = marray_cut(parIDs, :,:);
%     mcmc_3D(mcut0(:,:)', parnames(parIDs), '3 params');
%
%
%     mcmc_plot(marray_cut, parnames(:),...
%     'savematlabfig', figsave, 'savejpeg', false,...
%     'projdir', projdir, 'tstamp', tsToSave, 'extrafignamestring', 'BurnedIn');
%     mcut = marray_cut(parIDs, 1:end,ceil(end/4):end);
%     mcmc_3D(mcut(:,:)', parnames(parIDs), '3 params')
%



% We can literally pick any point in this cartesian product.
%     parRanges(parIDs, :) = [-1.16, -0.92 ;
%         -4 -0.5;
%         -15 -8.5]


%     parRanges([1, 2, 3], :) = [15.59 15.73 ;
%         -0.2966 0.1031;
%         8.385 8.622];
%     marray_cut = mcmc_cut(marray, [6,7,9], flipud((parRanges([6 7 9], :))'));
%
%     mcmc_plot(marray_cut(:, 1:end,ceil(end/4):end), parnames(:),...
%     'savematlabfig', figsave, 'savejpeg', false,...
%     'projdir', projdir, 'tstamp', tsToSave, 'extrafignamestring', 'BurnedIn');
%     mcut = marray_cut([6 7 9], 1:end,ceil(end/4):end);
%     mcmc_3D(mcut(:,:)', parnames([6 7 9]), 'RNA deg covariation')

% Here, we try to plot the trajectories for different sets of parameters.
% I will basically cut the parameters to different subsets, and plot the
% results to check if the lasR model is even working. if it is, then it
% will be clear that it just needs more simulation / fixing of the F rates
% and an estimation of the Kds.
%



% % % %
% % % %     paramIndices = (1:12);
% % % %     parRanges(paramIndices, :) = [...
% % % %         10.5 20.5; %1 pol_{Kd,tet}
% % % %         7 13; %2 pol_{Kd,lac}
% % % %         4 6.1; %3 pol % 5.9 in mcmc_info_ZSIFFL_training_full (from mtet phase 1, 2 etc)
% % % %         5.7 12; %4 Ribo % 5.9 in mcmc_info_ZSIFFL_training_full (from mtet phase 1, 2 etc)
% % % %        -20 20 ; %5 3OC12_{Kd} % correlated with plas_{tf, F}
% % % %         0 2; %6 3OC12_{F} % just gaussian.
% % % %         -20 20; %7 pol_{Kd,las}
% % % %         -2 2.2; %8 pol_{F,las}
% % % %         -20.5 20;%9 plas_{tf, Kd} <
% % % %         -20 20; % 10 -- plas-pol_{tf, Kd} <
% % % %         -2 5;% plas-pol_{tf, F}
% % % %         -2.5 1.5;% plas_{tf, F}
% % % %         ];
% % % %     marray_cut = mcmc_cut(marray, paramIndices, flipud((parRanges(paramIndices, :))'));
% % % %     size(marray_cut)

% actually i think it is better to just set te parameters using gaussian balls,
% and fix as many as possible, and jsut explore the Plas and OC12 Kds.
%
%     % mean and sd parameter.
%     gaussianMeanSD = ...
%         [16.57 0.1; ...%
%         10.5 0.1; ...% pol_{Kd,lac}
%         5.9 0.1; ...% pol
%         5.9 0.1; ...% Ribo
%         12 1; ...% 3OC12_{Kd}
%         0 1 ; ...% 3OC12_{F}
%         35 4; ...% pol_{Kd,las}
%         3 2; ...% pol_{F,las}
%         5 5; ...% plas_{tf, Kd} <
%         8 5; ...% plas-pol_{tf, Kd} <
%         0 2; ...% plas-pol_{tf, F}
%         0 4; ...%plas_{tf, F}
%         ];
%
%    marray_gauss =  cell2mat(arrayfun(@(mu, sig) mu + sig*randn(1, 100),...
%        gaussianMeanSD(:, 1), gaussianMeanSD(:, 2), ...
%        'UniformOutput', false));
%
%


% Plot trajectories.

% rebuild the master vector array, either via mcmc_cut or just using
% all estimated points.

%
%         mvarray = masterVecArray(marray_gauss, mai);
%          mvarray = masterVecArray(marray_cut, mai);
mvarray = masterVecArray(marray, mai);
%     clear marray
%
for miID = 1:length(mi)%1:
    currmi = mi(miID);
    dvStr = arrayfun(@num2str, currmi.dosedVals, 'UniformOutput', false);
    
    % titles array: # dose combs x (exp, sim) ==2 x # measured species.
    nms = length(currmi.measuredSpeciesIndex);
    ndc = size(currmi.dosedVals,2); %# dose combinations
    titls_array = cell(ndc, 2, nms);
    currdi = data_info(currmi.dataToMapTo);
    dn = currdi.dosedNames; %usually the last one will be the thing that is changing... 
    % can generalize to use all dose info if it becomes needed.
    for msID = 1:nms
        % for each measured species (ms), plot the trajectories over
        % all the doses, for both experiment and simulation
        if iscell(data_info(currmi.dataToMapTo).measuredNames{msID})
            if ischar( data_info(currmi.dataToMapTo).measuredNames{msID}{1})
                ms = data_info(currmi.dataToMapTo).measuredNames{msID}{1};
            else
                error('what is your data type?')
            end
        elseif ischar( data_info(currmi.dataToMapTo).measuredNames{msID})
            ms = data_info(currmi.dataToMapTo).measuredNames{msID};
        else
            error('what is your data type?')
        end
        for dcID = 1:ndc
            dosestr = [];
            
            for dnID = 1:length(dn)
                dosestr = ...
                    [dosestr ' ' dn{dnID} ' '...
                    num2str(currdi.dosedVals(dnID, dcID))];
            end
            % dose combination by (exp, sim) by measured species ID
            titls_array{dcID, 1, msID} = ['Exp ' ms dosestr];
            titls_array{dcID, 2, msID} = ['Sim ' ms dosestr];
        end
    end
    titls_array
    samplePoints = ceil(size(mvarray, 3) * [.95, 1]);
    
    marrayOrd = mvarray(currmi.paramMaps(currmi.orderingIx),:,samplePoints);
    
    if finafigmode == false
        mcmc_trajectories(currmi.emo, currdi, currmi, marrayOrd,...
            titls_array, {},...
            'SimMode', 'meanstd', 'separateExpSim', true,...
            'savematlabfig', figsave, 'savejpeg', jpgsave,...
            'projdir', projdir, 'tstamp', tsToSave,...
            'extrafignamestring', ...
            [extrastring num2str(miID) ' ' num2str(msID)]);%,'collateDoses', true,
    end
    
end
%
% end
marrayOrd(:,1:5,end)
[(1:42)' mvarray(:,1:5,end)]
% clear marrayOrd
% flagz = ones(26, 1);
% flagz([1 5 7 8 10 12 15 16 17 19 21 23 25 26]) = 0;
% [mvarray(:,1:6,end) log(cell2mat(activeNames2(:, 2))) flagz]

%% now make the collated trajectory figure and the endpoints figure for the paper.
% also, replot the markov chain trace plots and the posterior distribution
% plots here with proper font sizes.

% resimulate, but with just the matured degfp.
ms = {{'protein deGFP*'}};
for miID = 1:5
    currmi = mi(miID);
    currdi = data_info(currmi.dataToMapTo);
    tv = currdi.timeVector;
    marrayOrd = mvarray(currmi.paramMaps(currmi.orderingIx),:,samplePoints);
    dose  = currmi.dosedVals'
    vi = get(currmi.emo, 'ValueInfo')
    for i = 1:length(vi)
        vi(i)
    end
    currmi.dosedNames
    currmi.dosedVals
    
    dose
    
    [da{miID}, idxnotused{miID}] = ...
        simulatecurves(currmi.emo,marrayOrd(:,:)', 50, dose, tv, ms);
end


%% here


ylims = {[7.5, 2, 2, 5, 7.5];
    [7.5, 2, 2, 5, 7.5];
    [15 5 5 10 15]};


lengthToPlotArray = [24, 31, 61]; % 180 min = 24
for outercount = 1:length(lengthToPlotArray)
    lengthToPlot = lengthToPlotArray(outercount)
    
    figure
    set(0,'Units','normalized')
    figure
    set(gcf,'Units', 'normalized')
    set(gcf, 'Position', [0.05, 0.1, 0.9, 0.8])
    miToUse = [1 2 3 4 5];
    
    titleArray = {'pTet-deGFP DNA';
        'pLac-TetR DNA';
        'aTc';
        'pLac-deGFP DNA';
        '3OC12HSL'};
    
    legends = [{'4nM', '2nM', '1nM', '0.5nM', '250pM', '125pM', '0nM'}
        {'2nM', '.2nM', '.02nM', '2pM', '0.2pM', '20fM', '0nM'};
        {'10uM', '1uM', '0.1uM', '10nM', '1nM', '0.1nM', '0nM'}
        {'2nM', '1nM', '.5nM', '.25nM', '.125nM', '62.5pM', '0nM'};
        {'10uM', '1uM', '0.1uM', '10nM', '1nM', '0.1nM', '0nM'}];
    
    
    
    %     clear marray
    colorz = flipud(parula(length(legends(1, :))+2));
    for count = 1:length(miToUse)  %1:length(mi)%1:
        currmi = mi(miToUse(count));
        currdi = data_info(currmi.dataToMapTo);
        tv = currdi.timeVector;
        
        
        
        
        subplot(5, 2, (count-1)*2+1)
        
        for dID = 1:length(legends(count, :))
            
            plot(tv(1:lengthToPlot)/60,...
                0.001*mean(currdi.dataArray(1:lengthToPlot, 1, 1, dID),3),...
                'LineWidth', 1.5,...
                'Color', colorz(dID+2, :))
            
            hold on
        end
        
        title([titleArray{count} ' (Exp)'])
        ylabel('deGFP, uM', 'FontSize', 14)
        axis([0 (lengthToPlot-1)*8, 0 ylims{outercount}(count)])
        if count ==length(miToUse)
            xlabel('time, minutes', 'FontSize', 14)
        end
        ax = gca;
        ax.FontSize = 14;
        
        
        
        subplot(5, 2, (count-1)*2+2)
        for dID = 1:length(legends(count, :))
            
            plot(tv(1:lengthToPlot)/60,...
                0.001*mean(da{count}(1:lengthToPlot, 1, :, dID), 3),...
                'LineWidth', 1.5,...
                'Color', colorz(dID+2, :))
            
            hold on
        end
        
        title([titleArray{count} ' (Fit)'])
        legend(legends(count, :), 'Location', 'NorthWest', 'FontSize', 10)
        legend('boxoff')
        
        axis([0 (lengthToPlot-1)*8 0 ylims{outercount}(count)])
        if count ==length(miToUse)
            xlabel('time, minutes', 'FontSize', 14)
        end
        ax = gca;
        ax.FontSize = 14;
    end
    
    print([saveFinalFigs tsToSave '_fitting_traj_' num2str((lengthToPlot-1)*8) 'min'],'-depsc')
    print([saveFinalFigs tsToSave '_fitting_traj_' num2str((lengthToPlot-1)*8) 'min'],'-djpeg')
    
    print([saveFinalFigs tsToSave '_fitting_traj_' num2str((lengthToPlot-1)*8) 'min'],'-dpng')
    % create the endpoint plots
    % The first location is empty, for the IFFL schematic.
    figure
    exp_endpoints = zeros(size(legends)); % number of experiments by number of doses.
    exp_endpointsSD = zeros(size(legends));
    sim_endpoints = zeros(size(legends));
    sim_endpointsSD = zeros(size(legends));
    doseArray = zeros(size(legends));
    for i = 1:5
        doseArray(i, :) = data_info(i).dosedVals;
    end
    % doseArray(1, end) = 0.1/100; % 30c12
    % doseArray(2, end) = 0.0625/4; % lasR dna
    % doseArray(3, end) = 0.1/100; % aTc
    % doseArray(4, end) = 0.00001/100; % tetR dna
    % doseArray(5, end) = 0.125/4; % deGFP dna
    
    
    % build the data array
    for count = 1:length(miToUse)
        currmi = mi(miToUse(count));
        currdi = data_info(currmi.dataToMapTo);
        
        % get the experimental data endpoints
        for dID = 1:length(legends(count, :))
            exp_endpoints(count, dID) = ...
                0.001*mean(currdi.dataArray(lengthToPlot, 1, :, dID),3);
            exp_endpointsSD(count, dID) = ...
                std(0.001*currdi.dataArray(lengthToPlot, 1, :, dID),0,3);
            sim_endpoints(count, dID) = ...
                0.001*mean(da{count}(lengthToPlot, 1, :, dID), 3, 'omitnan');
            sim_endpointsSD(count, dID) = ...
                std(0.001*da{count}(lengthToPlot, 1, :, dID),0,  3, 'omitnan');
        end
    end
    %
    
    % Now generate the plots
    % close all
    
    %
    set(0,'Units','normalized')
    figure
    set(gcf,'Units', 'normalized')
    set(gcf, 'Position', [0.05, 0.1, 0.9, 0.8])
    
    
    xTickLabels = [{'4', '1', '2^{-2}',  '2^{-4}'};
        {'2', '.02',  '.0002',  '0'};
        {'10', '10^{-1}', '10^{-3}', '0'};
        {'2', '0.5', '0.125', '0'};
        {'10',  '10^{-1}',  '10^{-3}',  '0'}];
    titleArray = {'pTet deGFP DNA';
        'tetR DNA';
        'aTc';
        'pLac deGFP DNA';
        '3OC12'};
    xLab = {'pTet deGFP DNA';
        'tetR DNA';
        'aTc';
        'pLac deGFP DNA';
        '3OC12'};
    
    for count = 1:length(miToUse)
        subplot(2, 3, (count+1))
        ax = gca
        ax.XScale = 'log';
        errorbar(doseArray(count, :),...
            exp_endpoints(count, :),...
            exp_endpointsSD(count, :),...
            'LineWidth', 2);
        
        hold on
        
        errorbar(doseArray(count, :), ...
            sim_endpoints(count, :), ...
            sim_endpointsSD(count, :), ...
            'LineWidth', 2)
        ax = gca
        ax.XScale = 'log';
        hold on
        ax = gca
        ax.XTick = fliplr(doseArray(count, 1:2:end));
        ax.XTickLabel = fliplr(xTickLabels(count, :));
        ax.FontSize = 16
        
        %     title(titleArray{count})
        ylabel('deGFP, uM');
        xlabel(xLab{count});
        
    end
    
    print([saveFinalFigs tsToSave '_fitting_end_' num2str((lengthToPlot-1)*8) 'min'],'-depsc')
    print([saveFinalFigs tsToSave '_fitting_end_' num2str((lengthToPlot-1)*8) 'min'],'-djpeg')
    print([saveFinalFigs tsToSave '_fitting_end_' num2str((lengthToPlot-1)*8) 'min'],'-dpng')
    
end
%
