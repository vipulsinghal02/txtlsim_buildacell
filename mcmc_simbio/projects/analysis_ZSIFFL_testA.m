
% analysis_ZSIFFL_testA.m
% In this file, we simulate the IFFL for the full set of test experiments
% (5 in all). We will be drawing upon parameters from the parameter set of
% training_fullE, which is in itself a culmination of a lot of previous
% training. 
% We will define a new model, the IFFL, and an associated MCMC info, from
% which dosing info etc will be taken for mcmc_trajectories to plot the
% information. 

clear all
close all

% Set working directory to the txtlsim toolbox directory.
projdir = [pwd '/mcmc_simbio/projects/proj_ZSIFFL_testA'];
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
mcmc_info = mcmc_info_ZSIFFL_training_fullE(mtet, mlac, mlas);
di = data_ZSIFFL;
mi = mcmc_info.model_info;
ri = mcmc_info.runsim_info;
mai = mcmc_info.master_info;
% plot data from existing simulations.
tsIDtouse = 1;
plotflag = true;
switch tsIDtouse
    case 1 % this is after including the termination parameters.
        ts1 = '20190512_033129_1_1476';
        ts2 = '20190512_064547_1_1476';
        ts3 = '20190512_094712_1_1845';
        tstamp = {ts1 ts2};
        nIterID = {1:3 1:8 1:2};
        load([projdir '/simdata_' ts1 '/full_variable_set_' ts1 '.mat'], ...
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
    mcmc_plot(marray(:, 1:end,1:4:end), parnames(:),...
        'savematlabfig', figsave, 'savejpeg', jpgsave,...
        'projdir', projdir, 'tstamp', tsToSave,...
        'extrafignamestring', 'Without_transient');
% % % % %% 
%
        mcmc_plot(marray(:, 1:end,(end - 20):4:end), parnames(:),...
        'savematlabfig', figsave, 'savejpeg', jpgsave,...
        'projdir', projdir, 'tstamp', tsToSave, 'extrafignamestring', 'Burned_in');
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
