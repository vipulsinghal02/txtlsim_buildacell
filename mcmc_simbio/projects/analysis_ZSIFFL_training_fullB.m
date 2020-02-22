

% clear all
close all

% Set working directory to the txtlsim toolbox directory.
projdir = [pwd '/mcmc_simbio/projects/proj_ZSIFFL_trainingB'];
addpath(projdir)
sls = regexp(projdir, '/', 'split');
extrastring = sls{end};
jpgsave = true;
figsave = false;

% Load model, mcmc_info, and data_info.
% construct simbiology model object(s)
mtet = model_txtl_ptetdeGFP_pLactetR_aTc;
mlas = model_txtl_pLacLasR_pLasdeGFP;
mlac = model_txtl_pLacdeGFP;
% setup the mcmc_info struct
mcmc_info = mcmc_info_ZSIFFL_training_fullB(mtet, mlac, mlas);
di = data_ZSIFFL;
mi = mcmc_info.model_info;
ri = mcmc_info.runsim_info;
mai = mcmc_info.master_info;
% plot data from existing simulations.
tsIDtouse = 1;
plotflag = true;
switch tsIDtouse
    case 1 % this is after including the termination parameters. 
        ts1 = '20190508_001705_1_1845';
        ts2 = '20190508_213950_1_1107';
        ts3 = '20190509_024244_1_738';
        ts4 = '20190509_062932_1_738';
        tstamp = {ts1 ts2 ts3 ts4};
        nIterID = {1:7 1:9 1:6 1:2};
        load([projdir '/simdata_' ts1 '/full_variable_set_' ts1 '.mat'], ...
            'mi',...
            'mcmc_info', 'data_info',  'ri'); 
end
tsToSave = ts4;
mai.masterVector
marray_full = mcmc_get_walkers(tstamp,nIterID, projdir);
marray = marray_full(:,:,1:end);
clear marray_full
parnames = ...
    [{'pol_{term}'}
    {'Ribo_{term}'}
    {'pol'}
    {'Ribo'}
    {'3OC12_{Kd}'}
    {'pol_{Kd,las}'}
    {'plas_{tf, Kd}'}
    {'plas-pol_{tf, Kd}'}    ];
% activeNames(estParamsIX,:)
% 
% ans =
% 
%   12×3 cell array
% 
%     {'TXTL_PTET_RNAPbound_Kd'  }    {[1.5713e+07]}    {1×2 double}
%     {'TXTL_PLAC_RNAPbound_Kd'  }    {[3.6316e+04]}    {1×2 double}
%     {'RNAP'                    }    {[  365.0375]}    {1×2 double}
%     {'Ribo'                    }    {[  365.0375]}    {1×2 double}
%     {'TXTL_INDUCER_LASR_AHL_Kd'}    {[    0.1353]}    {1×2 double}
%     {'TXTL_INDUCER_LASR_AHL_F' }    {[    3.6693]}    {1×2 double}
%     {'TXTL_PLAS_RNAPbound_Kd'  }    {[    0.1353]}    {1×2 double}
%     {'TXTL_PLAS_RNAPbound_F'   }    {[    3.6693]}    {1×2 double}
%     {'TXTL_PLAS_TFBIND_Kd'     }    {[    0.1353]}    {1×2 double}
%     {'TXTL_PLAS_TFRNAPbound_Kd'}    {[    7.3891]}    {1×2 double}
%     {'TXTL_PLAS_TFRNAPbound_F' }    {[    3.6693]}    {1×2 double}
%     {'TXTL_PLAS_TFBIND_F'      }    {[    3.6693]}    {1×2 double}
%
% if plotflag
% % % % 
% % % % close all
% % % % % %    close all
% % % % %%     % Plot trace and corner (posterior distribution) plots
% % % %     mcmc_plot(marray(:, 1:end,(end-120):5:end), parnames(:),...
% % % %         'savematlabfig', figsave, 'savejpeg', jpgsave,...
% % % %         'projdir', projdir, 'tstamp', tsToSave,...
% % % %         'extrafignamestring', 'AllWalkers');
% % % % %% 
        mcmc_plot(marray(:, 1:end,1:end), parnames(:),...
        'savematlabfig', figsave, 'savejpeg', jpgsave,...
        'projdir', projdir, 'tstamp', tsToSave, 'extrafignamestring', 'Burned_in');
% % % % %

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
    
    
%         mvarray = masterVecArray(marray_gauss, mai);
%         mvarray = masterVecArray(marray_cut', mai);
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
[(1:42)' mvarray(:,1:5,end)]
clear marrayOrd
% flagz = ones(26, 1);
% flagz([1 5 7 8 10 12 15 16 17 19 21 23 25 26]) = 0;
% [mvarray(:,1:6,end) log(cell2mat(activeNames2(:, 2))) flagz]
