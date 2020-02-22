

% clear all
% clear all

%%
% Set working directory to the txtlsim toolbox directory.
projdir = [pwd '/mcmc_simbio/projects/proj_ZSIFFL_mtet'];

sls = regexp(projdir, '/', 'split');
extrastring = sls{end};
addpath(projdir)

jpgsave = true;
figsave = true;

% Load model, mcmc_info, and data_info.
mtet = model_txtl_ptetdeGFP_pLactetR_aTc;
mcmc_info = mcmc_info_ZSIFFL_mtet_phase1(mtet);
di = data_ZSIFFL;

mi = mcmc_info.model_info;
ri = mcmc_info.runsim_info;
mai = mcmc_info.master_info;

% plot data from existing simulations.
tsIDtouse = 3;
plotflag = true;
switch tsIDtouse
    case 1
        ts1 = '20190420_045130_1_92790';
        ts2 = '20190420_070950_1_30930';
        tstamp = {ts1 ts2};
        nIterID = {1:2, 1};
        load([projdir '/simdata_' ts1 '/full_variable_set_' ts1 '.mat'], ...
            'mi',...
            'mcmc_info', 'data_info', 'mai', 'ri');
    case 2
        ts1 = '20190420_045130_1_92790';
        ts2 = '20190420_071438_1_30930';
        ts3 = '20190420_071438_2_15465';
        ts4 = '20190420_155505_1_7732';
        ts5 = '20190420_155505_2_3093';
        ts6 = '20190420_155505_3_1546';
        ts7 = '20190420_155505_4_773';
        ts8 = '20190420_155505_5_309';
        ts9 = '20190420_155505_6_155';
        ts10 = '20190420_155505_7_77';
        ts11 = '20190420_155505_8_31';
        %
        %         tstamp = {ts1 ts2 ts3 ts4 ts5 ts6 ts7 ts8 ts9 ts10 ts11};
        %         nIterID = {1:2, 1:10 1:8 1:5 1:5 1:5 1:5 1:5, 1:5 1:5 1:5};
        
        tstamp = {ts4 ts5 ts6 ts7 ts8 ts9 ts10 ts11};
        nIterID = {1:5 1:5 1:5 1:5 1:5, 1:5 1:5 1:5};
        load([projdir '/simdata_' ts1 '/full_variable_set_' ts1 '.mat'], ...
            'mi',...
            'mcmc_info', 'data_info', 'mai', 'ri');
        
        %         conclusions from case 2:
        %     [{'pol_{Kd, tet}'     } OK
        %     {'rep_{Kd}'} expand lower bound from ~3 to maybe -5?
        %     {'rep_F' }  OK
        %     {'ATC_{Kd}'  } expand upper bound from 18 to 25?
        %     {'ATC_F'   } OK
        %     {'dim_{Kd}'        } expand lower bound from 7 to -3
        %     {'dim_F'         }
        
    case 3
        ts1 = '20190421_155749_1_773';
        ts2 = '20190421_155749_2_309';
        ts3 = '20190422_142534_1_773';
        ts4 = '20190422_214228_1_773';
        ts5 = '20190424_104801_1_773';
        
        tstamp = {ts1 ts2 ts3 ts4 ts5};
        nIterID = {1:10 1:2 1:4 1:7 1:14};
        load([projdir '/simdata_' ts1 '/full_variable_set_' ts1 '.mat'], ...
            'mi',...
            'mcmc_info', 'data_info', 'mai', 'ri');
end
tsToSave = ts5;
mai.masterVector


marray_full = mcmc_get_walkers(tstamp,nIterID, projdir);
marray = marray_full(:,:,1:end);
clear marray_full


parnames = ...
    [    {'rep_{Kd}'}
    {'ATC_{Kd}'  }
    {'dim_{Kd}'        }    ];
%
%%
if plotflag
    %     close all
    % Plot trace and corner (posterior distribution) plots
%     mcmc_plot(marray(:, 1:end,(end - 70):end), parnames(:),...
%         'savematlabfig', figsave, 'savejpeg', jpgsave,...
%         'projdir', projdir, 'tstamp', tsToSave, 'extrafignamestring', 'AllWalkers');
    %
    
    % {'rep_{Kd}'} -1.16, -0.92
    %
    
    parIDs = [1 2 3];
    parRanges(parIDs, :) = [-1.16, -0.92 ;
        -4 -0.5;
        -15 -8.5]
    
    marray_cut = mcmc_cut(marray(:, 1:end,(end - 70):end), parIDs, flipud((parRanges(parIDs, :))'));
%     
%     mcmc_plot(marray_cut(:, 1:end,ceil(end/4):end), parnames(:),...
%         'savematlabfig', figsave, 'savejpeg', false,...
%         'projdir', projdir, 'tstamp', tsToSave, 'extrafignamestring', 'BurnedIn');
    mcut = marray_cut(parIDs, 1:end,ceil(end/4):end);
%     mcmc_3D(mcut(:,:)', parnames(parIDs), '3 params')
    
    
    
    
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
    %
    %     paramIndices = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13];
    %     parRanges(paramIndices, :) = [...
    %         2.3 2.95; %1
    %         8.6 8.95; %2
    %         -10.1 -9.7; %3
    %         13 15; %4
    %         -1 5 ; %5
    %         15.59 15.73; %6
    %         -0.2966 0.1031; %7
    %         0.7 2.2; %8 pol
    %         8.385 8.622;%9 RNase
    %         3.18 3.69; % 10 -- TLcat
    %         -3 13.5;%RiboKd
    %         2 3;% Ribo term
    %         3 4.6 ];% Ribo
    %     marray_cut = mcmc_cut(marray, paramIndices, flipud((parRanges(paramIndices, :))'));
    %
    %     mcmc_plot(marray_cut(:, 1:end,ceil(end/4):end), parnames(:),...
    %     'savematlabfig', figsave, 'savejpeg', false,...
    %     'projdir', projdir, 'tstamp', tsToSave, 'extrafignamestring', 'BurnedIn');
    %
    %  CandidateParams = marray_cut(:,1:100:end,end)
    %     figure
    %     [C,lags,ESS]=eacorr(marray(:, :,1:end));%10000:end
    %     plot(lags,C,'.-',lags([1 end]),[0 0],'k');
    %     grid on
    %     xlabel('lags')
    %     ylabel('autocorrelation');
    %     text(lags(end),0,sprintf('Effective Sample Size (ESS): %.0f_ ',...
    %         ceil(mean(ESS))),'verticalalignment','bottom','horizontalalignment','right')
    %     title('Markov Chain Auto Correlation')
    %
    %% Plot trajectories.
    
    % rebuild the master vector array, either via mcmc_cut or just using
    % all estimated points.
    
    %     mvarray = masterVecArray(marray_cut, mai);
    mvarray = masterVecArray(marray, mai);
    clear marray
    
    for miID = 1:length(mi)
        currmi = mi(miID);
        samplePoints = ceil(size(mvarray, 3) * [.9, 1]);
        marrayOrd = mvarray(currmi.paramMaps(currmi.orderingIx),:,samplePoints);
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
                    dosestr = [dosestr ' ' dn{dnID} ' ' num2str(currdi.dosedVals(dnID, dcID))];
                end
                % dose combination by (exp, sim) by measured species ID
                titls_array{dcID, 1, msID} = ['Exp ' ms dosestr];
                titls_array{dcID, 2, msID} = ['Sim ' ms dosestr];
            end
        end
        titls_array
        samplePoints = ceil(size(mvarray, 3) * [.9, 1]);
        
        marrayOrd = mvarray(currmi.paramMaps(currmi.orderingIx),:,samplePoints);
        
        
        mcmc_trajectories(currmi.emo, currdi, currmi, marrayOrd,...
            titls_array, {},...
            'SimMode', 'meanstd', 'separateExpSim', true,...
            'savematlabfig', figsave, 'savejpeg', jpgsave,...
            'projdir', projdir, 'tstamp', tsToSave,...
            'extrafignamestring', [extrastring num2str(miID) ' ' num2str(msID)]);
    end
    
%     
%     %
%     marrayOrd = mvarray(mi(2).paramMaps(mi(2).orderingIx),:,samplePoints);
%     clear mvarray
%     titls = arrayfun(@num2str, mi(2).dosedVals, 'UniformOutput', false);
%     titls_array = cell(length(titls), 2, length(mi(2).measuredSpeciesIndex));
%     ms = {'MG aptamer', 'deGFP'};
%     for i = 1:length(mi(2).measuredSpeciesIndex)
%         for j = 1:length(titls)
%             titls_array(j, 1, i) = {[ms{i} ', ' titls{j} 'nM initial DNA, Exp data']};
%             titls_array(j, 2, i) = {[ms{i} ', ' titls{j} 'nM initial DNA, MCMC samples']};
%         end
%     end
%     mcmc_trajectories(mi(2).emo, data_info(mi(2).dataToMapTo), mi(2), marrayOrd,...
%         titls_array, {},...
%         'SimMode', 'meanstd', 'separateExpSim', true,...
%         'savematlabfig', figsave, 'savejpeg', jpgsave,...
%         'projdir', projdir, 'tstamp', tsToSave,...
%         'extrafignamestring', 'MGa_deGFP');
 end
% marrayOrd(:,1:5,end)
% clear marrayOrd
% flagz = ones(26, 1);
% flagz([1 5 7 8 10 12 15 16 17 19 21 23 25 26]) = 0;
% [mvarray(:,1:6,end) log(cell2mat(activeNames2(:, 2))) flagz]
