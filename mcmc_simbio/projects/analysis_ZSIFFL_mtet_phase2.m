

% clear all


% Set working directory to the txtlsim toolbox directory.
projdir = [pwd '/mcmc_simbio/projects/proj_ZSIFFL_mtet_phase2'];
addpath(projdir)
sls = regexp(projdir, '/', 'split');
extrastring = sls{end};
jpgsave = true;
figsave = false;

% Load model, mcmc_info, and data_info.
% construct simbiology model object(s)
mtet = model_txtl_ptetdeGFP_pLactetR_aTc;
mlac = model_txtl_pLacdeGFP;
% setup the mcmc_info struct
mcmc_info = mcmc_info_ZSIFFL_mtet_phase2(mtet, mlac);

di = data_ZSIFFL;

mi = mcmc_info.model_info;
ri = mcmc_info.runsim_info;
mai = mcmc_info.master_info;

% plot data from existing simulations.
tsIDtouse = 2;
plotflag = true;
switch tsIDtouse
    
    case 1
        ts1 = '20190427_170334_1_2058';
        ts2 = '20190427_170334_2_1029';
        
        tstamp = {ts1 ts2};
        nIterID = {1:10 1:3};
        load([projdir '/simdata_' ts1 '/full_variable_set_' ts1 '.mat'], ...
            'mi',...
            'mcmc_info', 'data_info', 'mai', 'ri'); 
        
    case 2
        ts1 = '20190428_142033_1_2058';
        ts2 = '20190428_142033_2_1029';
        ts3 = '20190429_083138_1_1029';
        ts4 = '20190429_200219_1_1029';
        ts5 = '20190429_200219_2_412';
        ts6 = '20190430_141254_1_412';
        ts7 = '20190430_141254_2_206';
        ts8 = '20190501_042829_1_617';
        ts9 = '20190501_123015_1_617';
        ts10 = '20190502_105714_1_412';
        ts11 = '20190503_073640_1_412';
        ts12 = '20190503_113414_1_412';
        ts13 = '20190503_170440_1_412';
        ts14 = '20190504_155258_1_412';
        tstamp = {ts1 ts2 ts3 ts4 ts5 ts6 ts7 ts8 ts9 ts10 ts11 ts12 ts13 ts14};
        nIterID = {1:10 1:2 1:4 1:5 1:3 1:5 1:5 1:4 1:8 1:11 1 1 11 13};
        load([projdir '/simdata_' ts1 '/full_variable_set_' ts1 '.mat'], ...
            'mi',...
            'mcmc_info', 'data_info', 'mai', 'ri'); 
end
tsToSave = ts14;
mai.masterVector

marray_full = mcmc_get_walkers(tstamp,nIterID, projdir);
marray = marray_full(:,:,1:end);
clear marray_full

    
parnames = ...
    [...
    {'tx_{cat}'}
    {'tl_cat'}
    {'tau'}
    {'pol_{Kd,tet}'}
    {'rep_{Kd}'}
    {'ATC_{Kd}'  }
    {'pol_{Kd,lac}'        }
    {'pol'}
    {'RNase'}
    {'Ribo'}];
%     {'TL_elong_glob'             }
%     {'AGTPdeg_time'              }
%     {'TXTL_PTET_RNAPbound_Kd'    }
%     {'TXTL_PTET_sequestration_Kd'}
%     {'TXTL_INDUCER_TETR_ATC_Kd'  }
%     {'TXTL_PLAC_RNAPbound_Kd'    }
%     {'RNAP'                      }
%     {'RNase'                     }
%     {'Ribo'                      }

% if plotflag
    mcmc_plot(marray(:, 1:end,(end-(100)):10:end), parnames(:),...
        'savematlabfig', figsave, 'savejpeg', jpgsave,...
        'projdir', projdir, 'tstamp', tsToSave, 'extrafignamestring', 'AllW_last150_thin10');
    %%
%         mcmc_plot(marray(:, 1:5:end,1:end), parnames(:),...
%         'savematlabfig', figsave, 'savejpeg', jpgsave,...
%         'projdir', projdir, 'tstamp', tsToSave, 'extrafignamestring', 'Burned_in');
%%

% {'rep_{Kd}'} -1.16, -0.92
% 

% parIDs = [1 2 3];
%     parRanges(parIDs, :) = [-1.16, -0.92 ; 
%         -4 -0.5;
%         -15 -8.5]
%       
%     marray_cut = mcmc_cut(marray(:, 1:end,(end - 70):end), parIDs, flipud((parRanges(parIDs, :))'));
%     
%     mcmc_plot(marray_cut(:, 1:end,ceil(end/4):end), parnames(:),...
%     'savematlabfig', figsave, 'savejpeg', false,...
%     'projdir', projdir, 'tstamp', tsToSave, 'extrafignamestring', 'BurnedIn');
%     mcut = marray_cut(parIDs, 1:end,ceil(end/4):end);
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
%     clear marray
    %
    for miID = 1:length(mi)
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
                    dosestr = [dosestr ' ' dn{dnID} ' ' num2str(currdi.dosedVals(dnID, dcID))];
                end
                % dose combination by (exp, sim) by measured species ID
                titls_array{dcID, 1, msID} = ['Exp ' ms dosestr];
                titls_array{dcID, 2, msID} = ['Sim ' ms dosestr];
            end
        end
        titls_array
        samplePoints = ceil(size(mvarray, 3) * [.95, 1]);
        
        marrayOrd = mvarray(currmi.paramMaps(currmi.orderingIx),:,samplePoints);
        
        
        mcmc_trajectories(currmi.emo, currdi, currmi, marrayOrd,...
            titls_array, {},...
            'SimMode', 'meanstd', 'separateExpSim', true,...
            'savematlabfig', figsave, 'savejpeg', jpgsave,...
            'projdir', projdir, 'tstamp', tsToSave,...
            'extrafignamestring', [extrastring num2str(miID) ' ' num2str(msID)]);
    end
    %
% end
marrayOrd(:,1:5,end)
clear marrayOrd
% flagz = ones(26, 1);
% flagz([1 5 7 8 10 12 15 16 17 19 21 23 25 26]) = 0;
% [mvarray(:,1:6,end) log(cell2mat(activeNames2(:, 2))) flagz]
