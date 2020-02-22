


% In this file we try to find a set of parameter values that allow the
% model from model_dsg2014_regen to fit the data from data_dsg2014_full.
% Vipul Singhal, 2019

% simulate txtl model with custom parameter values, and look at the species
% plots as specified by mcmc_info object.

% Set working directory to the txtlsim toolbox directory.
projdir = [pwd '/mcmc_simbio/projects/proj_acs_dsg2014_regen_C1'];
addpath(projdir)

% Load model, mcmc_info, and data_info.
mobj = model_dsg2014_regen;
mcmc_info = mcmc_info_dsg2014_regen_C1(mobj);
di = data_dsg2014_full;

mi = mcmc_info.model_info;
ri = mcmc_info.runsim_info;
mai = mcmc_info.master_info;

% plot data from existing simulations.
tsIDtouse = 3;
plotflag = true;
switch tsIDtouse
    case 1 % these were all messed up because the master vector was a bunch of zeros (log(1)'s)
        % done on AWS, LHS initialization.
        ts1temp = '20190208_075812_1_217593'; % 1800 steps with ladder [1.4 1.2 1]*1.2
        ts2temp = '20190208_075812_2_108797';
        ts3temp = '20190208_075812_3_21759';
        
        nIterID = {1:7, 1:7, 1:6};
        tstamp = {ts1temp ts2temp ts3temp};
        load([projdir '/simdata_' ts1temp '/full_variable_set_' ts1temp '.mat'], ...
            'mi',...
            'mcmc_info', 'data_info', 'mai', 'ri');
    case 2 % feb 9, 3.10 pm. step 1.15, temp 0.05 (1088), 1200 walkers, 33 steps per iter, 20 iter. 
        ts1 = '20190209_070634_1_1088';
        nIterID = {1:13}
        tstamp = {ts1};
        load([projdir '/simdata_' ts1 '/full_variable_set_' ts1 '.mat'], ...
            'mi',...
            'mcmc_info', 'data_info', 'mai', 'ri');
    case 3
        ts1 = '20190209_091205_1_231';
        nIterID = {1:5}
        load([projdir '/simdata_' ts1 '/full_variable_set_' ts1 '.mat'], ...
            'mi',...
            'mcmc_info', 'data_info', 'mai', 'ri');
end

mai.masterVector


marray = mcmc_get_walkers(tstamp,nIterID, projdir);
parnames = ...
    [{'TX_{cat}'    }
    {'pol_{Kd}'     }
    {'pol_{term}'   }
    {'n_{Kd1}'      }
    {'n_{Kd2}'      }
    {'RNAse_{Kd}'   }
    {'pol'          }
    {'RNase'        }
    {'TL_{cat}'     }
    {'Ribo_{Kd}'    }
    {'aa_{Kd}'      }
    {'TL_n_{Kd}'    }
    {'Ribo_{term}'  }
    {'Ribo'         }];
%

if plotflag
    close all
    % Plot trace and corner (posterior distribution) plots
mcmc_plot(marray(:, 1:50:end,:), parnames(:));
%     mcmc_plot(marray(1:6, 1:20:end,1:end), parnames(1:6));
%     mcmc_plot(marray(6+(1:8), 1:20:end,1:end), parnames(6+(1:8)));

    figure
    [C,lags,ESS]=eacorr(marray(:, :,:));%10000:end
    plot(lags,C,'.-',lags([1 end]),[0 0],'k');
    grid on
    xlabel('lags')
    ylabel('autocorrelation');
    text(lags(end),0,sprintf('Effective Sample Size (ESS): %.0f_ ',...
        ceil(mean(ESS))),'verticalalignment','bottom','horizontalalignment','right')
    title('Markov Chain Auto Correlation')
    %
    titls = arrayfun(@num2str, mi(1).dosedVals, 'UniformOutput', false);
    titls_array = cell(length(titls), 1, length(mi(2).measuredSpeciesIndex));
    for i = 1:length(mi(2).measuredSpeciesIndex)
        for j = 1:length(titls)
            titls_array(j, 1, i) = titls(j);
        end
    end
    mvarray = masterVecArray(marray, mai);
    samplePoints = ceil(size(mvarray, 3) * [.8, 1]);
    %
    marrayOrd = mvarray(mi(1).paramMaps(mi(1).orderingIx),:,samplePoints);
    mcmc_trajectories(mi(1).emo, data_info(mi(1).dataToMapTo), mi(1), marrayOrd,...
        titls', {},...
        'SimMode', 'meanstd');%, 'savematlabfig', false, 'savejpeg', false
    
    marrayOrd = mvarray(mi(2).paramMaps(mi(2).orderingIx),:,samplePoints);
    mcmc_trajectories(mi(2).emo, data_info(mi(2).dataToMapTo), mi(2), marrayOrd,...
        titls_array, {},...
        'SimMode', 'meanstd');
end
marrayOrd(:,1:5,end)
flagz = ones(26, 1);
flagz([1 5 7 8 10 12 15 16 17 19 21 23 25 26]) = 0;
[mvarray(:,1:6,end) log(cell2mat(activeNames2(:, 2))) flagz]

% plot trajectories
lgds = {};
% aaa =...
%     [{[  2.1904]}    {[  2.1904]}    {'TX_elong_glob'     }
%     {[  6.1017]}    {[  6.1017]}    {'TL_elong_glob'     }
%     {[  6.7502]}    {[  6.7502]}    {'AGTPdeg_time'      }
%     {[ -8.9816]}    {[ -8.9816]}    {'AGTPdeg_rate'      }
%     {[ -1.9799]}    {[ -1.9799]}    {'TXTL_PROT_deGFP_?'}
%     {[ -1.4155]}    {[ -1.4155]}    {'TXTL_UTR_UTR1_Kd'  }
%     {[  1.0949]}    {[  1.0949]}    {'TXTL_UTR_UTR1_F'   }
%     {[  1.8325]}    {[  1.8325]}    {'TXTL_P70_RNAPbou?'}
%     {[  4.5064]}    {[  4.5064]}    {'TXTL_P70_RNAPbou?'}
%     {[  -2.297]}    {[  -2.297]}    {'TXTL_RNAPBOUND_T?'}
%     {[0.096026]}    {[0.096026]}    {'TXTL_NTP_RNAP_1_Kd'}
%     {[  2.1533]}    {[  2.1533]}    {'TXTL_NTP_RNAP_1_F' }
%     {[  2.2435]}    {[  2.2435]}    {'TXTL_NTP_RNAP_2_Kd'}
%     {[-0.83904]}    {[-0.83904]}    {'TXTL_NTP_RNAP_2_F' }
%     {[  3.4805]}    {[  3.4805]}    {'TL_AA_Kd'          }
%     {[-0.22471]}    {[-0.22471]}    {'TL_AA_F'           }
%     {[  5.6524]}    {[  5.6524]}    {'TL_AGTP_Kd'        }
%     {[ 0.97604]}    {[ 0.97604]}    {'TL_AGTP_F'         }
%     {[-0.95634]}    {[-0.95634]}    {'TXTL_RIBOBOUND_T?'}
%     {[  8.5748]}    {[  8.5748]}    {'TXTL_RNAdeg_Kd'    }
%     {[ -0.5005]}    {[ -0.5005]}    {'TXTL_RNAdeg_F'     }
%     {[ -5.1466]}    {[ -5.1466]}    {'TXTL_RNAdeg_kc'    }
%     {[  2.2309]}    {[  2.2309]}    {'RNAP'              }
%     {[  6.7048]}    {[  6.7048]}    {'Ribo'              }
%     {[  4.5145]}    {[  4.5145]}    {'RNase'             }];
%%
activeNames1 = {...
    'TXTL_RNAdeg_Kd'                    2000        [100 10000]
    'TXTL_RNAdeg_F'                     0.02        [0.01 100]
    'TXTL_RNAdeg_kc'                    0.0028      [1e-4 1]
    'RNase'                             100         [10 1000]};

simulatecurves(...
    mi(1).emo,...
    log(cell2mat(activeNames1(:,2))'),...
    1,...
    di(mi(1).dataToMapTo).dosedVals', ...
    di((mi(1).dataToMapTo)).timeVector, ...
    mi(1).measuredSpecies)

titls = arrayfun(@num2str, mi(1).dosedVals, 'UniformOutput', false);
mcmc_trajectories(mi(1).emo,...
    di(mi(1).dataToMapTo),...
    mi(1),...
    log(cell2mat(activeNames1(:,2))'),...
    titls',...
    {},...
    'SimMode', 'curves')
%%
figure
debugda = simulatecurves(...
    mi(1).emo,...
    log(cell2mat(activeNames1(:,2))'),...
    1,...
    di(mi(1).dataToMapTo).dosedVals', ...
    di((mi(1).dataToMapTo)).timeVector, ...
    mi(1).measuredSpecies);

for i = 1:size(debugda, 4)
    subplot(size(debugda, 4), 1, i)
    plot(di((mi(1).dataToMapTo)).timeVector, debugda(:,1,1,i));
    title(di(mi(1).dataToMapTo).dosedVals(i));
end
%%
debugsd = simulate(mi(1).emo);
speciesnames = debugsd.DataNames;
figure
for i = 1:30
    subplot(6, 5, i)
    currdata = simulatecurves(...
        mi(1).emo,...
        log(cell2mat(activeNames1(:,2))'),...
        1,...
        di(mi(1).dataToMapTo).dosedVals', ...
        di((mi(1).dataToMapTo)).timeVector, ...
        {speciesnames(i)});
    plot(di((mi(1).dataToMapTo)).timeVector, currdata(:,1,1,1));
    title(speciesnames{i})
    
    
end


%%
    activeNames2 = {... % changes made to ranges on feb 8, 2019. setting parameters based on 
        ...% posterior plots. 
        'TX_elong_glob'                     exp(3.4)      [0.5 300]        % 1
        'AGTPdeg_time'                      exp(9.57)   [1800 42000] % set to exp(9.57)
        'AGTPdeg_rate'                      0.0002      [1e-5 1e-2] % set from before
        'AGTPreg_ON'                        0.02        [0.005 0.2]        %4 % set from before
        'TXTL_P70_RNAPbound_Kd'             200         [0.1 1e6]  % 
        'TXTL_P70_RNAPbound_F'              exp(1.5)    [1e-5 300] % set to exp(1.5)
        'TXTL_RNAPBOUND_TERMINATION_RATE'   0.15        [1e-4 100]         % 7
        'TXTL_NTP_RNAP_1_Kd'                exp(12.5)      [1 1e6]
        'TXTL_NTP_RNAP_1_F'                 exp(0)      [1e-5 100] % set to 1
        'TXTL_NTP_RNAP_2_Kd'                exp(3)         [0.1 1e7]
        'TXTL_NTP_RNAP_2_F'                 exp(0)      [1e-6 1000]        %11 % set to 1
        'TXTL_RNAdeg_Kd'                    exp(7.6)        [100 1e5] 
        'TXTL_RNAdeg_F'                     1        [0.01 10000] % set to 1 (1 is right in the middle of the broad posterior density, and so not entirely arbitrary. )
        'TXTL_RNAdeg_kc'                    exp(-5.4)   [1e-4 1]  %set to exp(-5.4)
        'RNAP'                              exp(5.8)         [5 5000] % 15
        'RNase'                             exp(5.3)         [1 10000]
        'TL_elong_glob'                     20          [0.1 500]
        'TXTL_PROT_deGFP_MATURATION'        0.0023      [0.0002 0.02] %18 % set from before
        'TXTL_UTR_UTR1_Kd'                  exp(4.5)          [0.05 1e5]
        'TXTL_UTR_UTR1_F'                   exp(-0.2)   [1e-5 100] % set to exp(-0.2)
        'TL_AA_Kd'                          exp(7.5)      [.1 1e6] % 21
        'TL_AA_F'                           exp(-0.3)   [1e-5 20]   % set to exp(-0.3)
        'TL_AGTP_Kd'                        exp(7.5)     [.1 1e7] % 23
        'TL_AGTP_F'                         exp(-1.2)   [1e-5 100]  %set to exp(-1.2)
        'TXTL_RIBOBOUND_TERMINATION_RATE'   exp(5)          [0.1 20000]
        'Ribo'                              exp(3.5)         [10 10000]}; %26 

    % these are unordered parameters, from the fits from case 4, from
    % regen_A. (file: manual_txtlsim_params)
    close all
    
    pix = 34
%     mvarray load from case 4 of regen A analysis file:
%     manual_txtlsim_parameters.
 paramsToTest = mvarray(:,1:5:end,end);  
 
   

simulatecurves(...
    mi(1).emo,...
    log(cell2mat(activeNames2(mi(1).paramMaps,2))'),...
    1,...
    di(mi(1).dataToMapTo).dosedVals', ...
    di((mi(1).dataToMapTo)).timeVector, ...
    mi(1).measuredSpecies)


titls = arrayfun(@num2str, mi(1).dosedVals, 'UniformOutput', false);
mcmc_trajectories(mi(1).emo,...
    di(mi(1).dataToMapTo),...
    mi(1),...
    (paramsToTest(mi(1).paramMaps(mi(1).orderingIx,1),pix)'),...
    titls',...
    {},...
    'SimMode', 'curves')
%

titls = arrayfun(@num2str, mi(2).dosedVals, 'UniformOutput', false);
titls_array = cell(length(titls), 1, length(mi(2).measuredSpeciesIndex));
for i = 1:length(mi(2).measuredSpeciesIndex)
    for j = 1:length(titls)
        titls_array(j, 1, i) = titls(j);
    end
end
mcmc_trajectories(mi(2).emo,...
    di(mi(2).dataToMapTo),...
    mi(2),...
    (paramsToTest(mi(2).orderingIx,pix)'),...
    titls_array,...
    {},...
    'SimMode', 'curves')



%%
% sets that do well
close all
clc
activeNames2(:,1)
found1 = [       5.5354
         9.57
      -8.5172
       -3.912
        9.514
          1.5
       3.3005
       2.9459
            0
       13.997
            0
        9.237
            0
         -5.4
       2.4419
       6.4899
       5.1219
      -6.0748
       11.189
         -0.2
       6.5566
         -0.3
       14.509
         -1.2
        5.398
       2.7081]; % mrna expression shape is fantastic. expression level is a bit high. manually adjusting: 
   mod1 = [      ...
       4.9354 % tx elong, verified in mcmc_info_dsg2014_regen_D.m on Aug 26, 2019
         9.17 % AGTPdeg_time, verified in mcmc_info_dsg2014_regen_D.m on Aug 26, 2019
      -9.5172 % AGTPdeg_rate, verified in mcmc_info_dsg2014_regen_D.m on Aug 26, 2019
       -3.912 % AGTPreg_ON, verified in mcmc_info_dsg2014_regen_D.m on Aug 26, 2019
        9.514 % 'TXTL_P70_RNAPbound_Kd', verified in mcmc_info_dsg2014_regen_D.m on Aug 26, 2019
          1.5 % 'TXTL_P70_RNAPbound_F', verified in mcmc_info_dsg2014_regen_D.m on Aug 26, 2019
       3.3005 % rnap termination rate, verified in mcmc_info_dsg2014_regen_D.m on Aug 26, 2019
       2.9459%'TXTL_NTP_RNAP_1_Kd' , verified in mcmc_info_dsg2014_regen_D.m on Aug 26, 2019         
            0%'TXTL_NTP_RNAP_1_F' , verified in mcmc_info_dsg2014_regen_D.m on Aug 26, 2019     
       13.997%'TXTL_NTP_RNAP_2_Kd' , verified in mcmc_info_dsg2014_regen_D.m on Aug 26, 2019        
            0%'TXTL_NTP_RNAP_2_F' , verified in mcmc_info_dsg2014_regen_D.m on Aug 26, 2019     
        9.237%'TXTL_RNAdeg_Kd' , verified in mcmc_info_dsg2014_regen_D.m on Aug 26, 2019         %freeze
            0%'TXTL_RNAdeg_F' , verified in mcmc_info_dsg2014_regen_D.m on Aug 26, 2019      % freeze
         -4.4%'TXTL_RNAdeg_kc' , verified in mcmc_info_dsg2014_regen_D.m on Aug 26, 2019         %freeze
       1.4419%'RNAP' , verified in mcmc_info_dsg2014_regen_D.m on Aug 26, 2019       % rnap
       6.4899%'RNase' , verified in mcmc_info_dsg2014_regen_D.m on Aug 26, 2019      %freeze
       0.5219%'TL_elong_glob' , verified in mcmc_info_dsg2014_regen_D.m on Aug 26, 2019     
      -6.0748%'TXTL_PROT_deGFP_MATURATION' , verified in mcmc_info_dsg2014_regen_D.m on Aug 26, 2019        
       11.189%'TXTL_UTR_UTR1_Kd' , verified in mcmc_info_dsg2014_regen_D.m on Aug 26, 2019      
         -0.2%'TXTL_UTR_UTR1_F' , verified in mcmc_info_dsg2014_regen_D.m on Aug 26, 2019       
       6.5566%'TL_AA_Kd' , verified in mcmc_info_dsg2014_regen_D.m on Aug 26, 2019      
         -0.3%'TL_AA_F' , verified in mcmc_info_dsg2014_regen_D.m on Aug 26, 2019       
       14.509%'TL_AGTP_Kd' , verified in mcmc_info_dsg2014_regen_D.m on Aug 26, 2019        
         -1.2%'TL_AGTP_F' , verified in mcmc_info_dsg2014_regen_D.m on Aug 26, 2019     
        5.398%'TXTL_RIBOBOUND_TERMINATION_RATE' , verified in mcmc_info_dsg2014_regen_D.m on Aug 26, 2019        % ribo termination rate
       7.3081 ]; %'Ribo' , verified in mcmc_info_dsg2014_regen_D.m on Aug 26, 2019 
           
% wow there is a trade off between the ribo number and the tl elongation
% rate. thats awesome. there is a linear vs saturation regime argument that
% needs to be made here. 
% if protein expression at high dna fits, but at low dna the sim are much
% higher than the exp, then in the model the ribo is saturated, and the tl
% elongation rate is too high. if you increase the ribo count, all the
% expression levels go super high, but now the ribo is not saturated, and
% so the high dna protein conc goes much higher than the low dna protein
% conc. --> now lower the elongation rate, and everything scaled back down
% proportionately. this applied to both rna and protein. this along with
% non identifiability is a short note / section in this paper. --> the
% relative expression at diff dna conc allows for identifiability between
% the ribo count and the tl elongation rate, all else fixed. this, adding
% one extra trajectory from 1 -> 2 trajectories (ie 2 diff dna conc)
% increases the identifiability. 

% parameter tradeoffs and non identifiability paper -- can describe a
% general case of this above. for a tx, tl and rna deg example system. 

% 
%
%
%
% ok, so with these parameter values the rna deg fits almost perfectly. the
% transcription and translation fit well for the highest dose. So we can
% start using these values as initial points. 
% also some observations/guesses: NTP and aa binding F rates and Kds
% probaby do not matter too much, and just serve to make the model much
% harder to fit. just get rid of these by fixing them to the values in mod1
% above. 
% 
% so set up the following sims: 
% for all sims: 
 
% 
%     {'AGTPreg_ON'                     }
%     {'TXTL_P70_RNAPbound_F'           }
%     {'TXTL_NTP_RNAP_1_Kd'             }
%     {'TXTL_NTP_RNAP_1_F'              }
%     {'TXTL_NTP_RNAP_2_Kd'             }
%     {'TXTL_NTP_RNAP_2_F'              }
%     {'TXTL_RNAdeg_F'                  }    
%     {'TXTL_UTR_UTR1_F'                }
%     {'TL_AA_Kd'                       }
%     {'TL_AA_F'                        }  
%     {'TL_AGTP_Kd'                     }
%     {'TL_AGTP_F'                      }  
%     {'TXTL_PROT_deGFP_MATURATION'     }
% are fixed at the values above. 

% 1. The minimal sim (regen_D) 
% also fix:
%     {'TXTL_RNAdeg_Kd'                 }
% and the AGTP params:
%     {'AGTPdeg_time'                   }
%     {'AGTPdeg_rate'                   }
% leaving only 


%     {'TX_elong_glob'                  }
%     {'TXTL_P70_RNAPbound_Kd'          }
%     {'TXTL_RNAPBOUND_TERMINATION_RATE'}
%     {'TXTL_RNAdeg_kc'                 }
%     {'RNAP'                           }
%     {'RNase'                          }
%     {'TL_elong_glob'                  }
%     {'TXTL_UTR_UTR1_Kd'               }
%     {'TXTL_RIBOBOUND_TERMINATION_RATE'}
%     {'Ribo'                           }
% 10 parameters to estimate. 
% can do this for both acsdsg and vnprl (vnprl has only 30nM ish for rna
% production. everything else is about the same, maybe use the figure in
% that paper to also scale the protein production levels by the one hour
% measurement. 

% 2. regen_E:
% fix the AGTP params:
%     {'AGTPdeg_time'                   }
%     {'AGTPdeg_rate'                   }
% leaving only 

%     {'TX_elong_glob'                  }
%     {'TXTL_P70_RNAPbound_Kd'          }
%     {'TXTL_RNAPBOUND_TERMINATION_RATE'}
%     {'RNAP'                           }
%     {'TL_elong_glob'                  }
%     {'TXTL_UTR_UTR1_Kd'               }
%     {'TXTL_RIBOBOUND_TERMINATION_RATE'}
%     {'Ribo'                           }
%     {'TXTL_RNAdeg_Kd'                 }
%     {'TXTL_RNAdeg_kc'                 }
%     {'RNase'                          }
% 11 parameters to estimate. RNA deg here allows the algorithm to tune the RNA timescale.  

% 3. regen_F:
% estimate all 14 remaining params:
%     {'AGTPdeg_time'                   }
%     {'AGTPdeg_rate'                   }
%     {'TX_elong_glob'                  }
%     {'TXTL_P70_RNAPbound_Kd'          }
%     {'TXTL_RNAPBOUND_TERMINATION_RATE'}
%     {'RNAP'                           }
%     {'TL_elong_glob'                  }
%     {'TXTL_UTR_UTR1_Kd'               }
%     {'TXTL_RIBOBOUND_TERMINATION_RATE'}
%     {'Ribo'                           }
%     {'TXTL_RNAdeg_Kd'                 }
%     {'TXTL_RNAdeg_kc'                 }
%     {'RNase'                          }


% regen_G: (rna deg topology removed. 
%     {'TX_elong_glob'                  }
%     {'TXTL_P70_RNAPbound_Kd'          }
%     {'RNAP'                           }
%     {'TL_elong_glob'                  }
%     {'TXTL_UTR_UTR1_Kd'               }
%     {'Ribo'                           }

% Do all of these at c5.9xlarge or m5.24xlarge, with step 1.18, temperature 0.001, 1200
% walkers, and 20 iter at 100*1200 points per iter. 

% and even among these 6, prioritize regen_D on both the ACS dsg data and
% the vnprl data. BOSS -- just do those two to start, and email RMM the
% results! 
% Plot with separate sim and exp, with max normalized y axes
    
   titls = arrayfun(@num2str, mi(1).dosedVals, 'UniformOutput', false);
mcmc_trajectories(mi(1).emo,...
    di(mi(1).dataToMapTo),...
    mi(1),...
    (mod1(mi(1).paramMaps(mi(1).orderingIx,1))'),...
    titls',...
    {},...
    'SimMode', 'curves')
%

titls = arrayfun(@num2str, mi(2).dosedVals, 'UniformOutput', false);
titls_array = cell(length(titls), 1, length(mi(2).measuredSpeciesIndex));
for i = 1:length(mi(2).measuredSpeciesIndex)
    for j = 1:length(titls)
        titls_array(j, 1, i) = titls(j);
    end
end
mcmc_trajectories(mi(2).emo,...
    di(mi(2).dataToMapTo),...
    mi(2),...
    (mod1(mi(2).orderingIx)'),...
    titls_array,...
    {},...
    'SimMode', 'curves')
