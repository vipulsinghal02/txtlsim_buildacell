function [mcmc_info, varargout] = mcmc_info_ZSIFFL_training_fullB(mtet, mlac, mlas)
% !! todo -- set the parameters estimated from phase 1. then this will be
% ready to run. 
% Phase 2 of the estimation procedure of the tetR system parameters from
% Zach's tetR data.
% 
% The topologies involved in this estimataion probem are: 
% plac - UTR1 - tetR, ptet - UTR1 - deGFP, aTc; This has three
% topologies, the constitutive expression geometry, the repression
% geometry, and the aTc induction geometry.
%
% A lot of the core parameters are set from those estimated from the
% constitutive expression fits to the data in the VNPRL2011 and ACS 2014 papers.
% 
% In phase 1, we estimated 7 parameters: 
%     {'TXTL_PTET_RNAPbound_Kd'    }
%     {'TXTL_PTET_sequestration_Kd'}
%     {'TXTL_PTET_sequestration_F' }
%     {'TXTL_INDUCER_TETR_ATC_Kd'  }
%     {'TXTL_INDUCER_TETR_ATC_F'   }
%     {'TXTL_DIMER_tetR_Kd'        }
%     {'TXTL_DIMER_tetR_F'         }
%
% In this phase, we will fix the values of
%     {'TXTL_PTET_sequestration_F' }
%     {'TXTL_INDUCER_TETR_ATC_F'   }
%     {'TXTL_DIMER_tetR_Kd'        }
%     {'TXTL_DIMER_tetR_F'         }
% 
% and estimate:
% 
%     'TX_elong_glob'                      , exp(2.6),   [exp(0) exp(5)]          %1 % from est params above
%     'TL_elong_glob'                      , exp(3.5),   [exp(0) exp(6)]          %2% from est params above
%     'AGTPdeg_time'                       , exp(8.8),   [exp(6) exp(11)]         %3 % from est params above
%     'TXTL_PTET_RNAPbound_Kd'             , exp(14),   [exp(0) exp(17)]          %7 % TO BE ESTIMATED HERE
%     'TXTL_PTET_sequestration_Kd'         , exp(12),   [exp(3) exp(15)]          %11 TO BE ESTIMATED HERE
%     'TXTL_INDUCER_TETR_ATC_Kd'           , exp(13),   [exp(0) exp(18)]          %16 % TO BE ESTIMATED HERE
%     'TXTL_PLAC_RNAPbound_Kd'             , exp(13.8),   [exp(5) exp(17)]        %21 from est params above
%     'RNAP'                               , exp(1.4419),  [exp(-1) exp(4)]       %31 31% from est params above
%     'RNase'                              , exp(8.5),  [exp(5) exp(10)]          %32 from est params above
%     'Ribo'                               , exp(3.75),  [exp(1) exp(6)]          %33 % from est params above
%
% and we also include the pLac data into the estimation problem, since the
% plac part is involved. 
%
%
% mcmc_info has the following substructures:
%
% runsim_info:  information on the mcmc algorithm parameters
% model_info:   array of models, and associated properties like parameters,
%               and the matrices of indices from the master vector
%               to the model parameters.
% master_info:  contains the master vector, and a spec for which parameters
%               get estimated.
%
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


% Some human readable descriptive text.
ptetinfo = ['This is the ptet constitutive expression circuit. \n '];

tetrinfo = ['This is the tetR repressing ptet circuit. \n '];

atcinfo = ['This is the aTc inducing (derepressing) the tetR ptet repression circuit. \n '];

placinfo = ['This is the plac constitutive expression circuit. \n '];

plasinfo = ['This is the plas - 3OC12HSL induction circuit. \n '];
% activeNames has the mRNA parameters and the protein parameters.
% first half (up to RNase) are TX and the rest are TL.
% TX params are fixed from previous sims.
%
% ordering requirements:
% ensure that the following two orderings match up:
% activeNames(orderingIX) == masterVector(paramMaps(orderingIX))
%
% ie, activeNames == masterVector(paramMaps)
%
% This gets satisfied when two conditions hold:
%
% The fixed parameters in the master vector must be arranged to that
% for every paramMap and every corresponding activeNames list, the
% fixed params subset of the elements gets mapped correctly.
%
% for the estimaed parameters, again, the estimated parameters need to
% populate the master vector in a way such that the condition
% activeNames == masterVector(paramMaps) holds for all the activeNames
% arrays (each topology will have one), and for each paramMap column
% (geometry) for each topology.
%
%
% of course, the masterVector is built as follows:
% masterVector(estparamIX) == logp
% masterVector(fixedParams) == [marray(:, end-2); marray(:, end-1); marray(:, end);]
%
% So this is all a bit complicated...
% Basically we need to make sure that after we build master vector from
% the fixed parameters (from the previous simulations), when we access
% them using paramMaps, we get the ones corresponding to the names in
% activeNames.


% names of the parameters and species to set allow for setting in the
% exported model. These are both the set parameters and the to estimate
% parameters.
% Here for example we have for the mrna deg sim we only care about
% setting the rna deg parameters.
%%

% SEE FILE ZachIFFL_estimation_strategy.txt for the overall strategy. 

% first we set up the master vector
%
% then we will set the parameters to fix, and the remaining ones to
% estimate.
%
% Then we will define paramMaps that will be used to distribute the
% master vector among the different topologies.

% parameters of the master vector we fix and estimate

%     The following params are the ones estimaed from the last set of
%     restrictions in the file analysis_vnprl_F2.m
% parnames =
%
%   13×1 cell array
%
%     {'TX_{cat}'    } 1
%     {'\tau_{atp}'  } 3
%     {'\delta_{atp}'} 5
%     {'pol_{Kd}'    } 21 % plac Kd NOT PTET Kd, the Ptet KD gets estimated. 
%     {'pol_{term}'  } 23
%     {'RNAse_{Kd}'  } 15
%     {'RNAse_{cat}' } 30
%     {'pol'         } 31
%     {'RNase'       } 32
%     {'TL_{cat}'    } 2
%     {'Ribo_{Kd}'   } 6
%     {'Ribo_{term}' } 28
%     {'Ribo'        } 33


% The following were parameters esitmated in vnprl_F2, and were extracted
% carefully from the accompanying analysis script. see the code there to
% understand how they were picked. Need to document this before I forget. 
EstimatedParams =[...
    2.5234    2.4464    2.6231    2.5976    2.4650    2.4277    2.4806    2.6010    2.4301    2.5499    2.6991    2.4356
    8.8054    8.9024    8.8097    8.8014    8.8483    8.8642    8.8394    8.8427    8.8616    8.7634    8.8621    8.7910
  -10.0179   -9.7873   -9.9376  -10.0048   -9.9056   -9.8889   -9.9100   -9.9631   -9.9827  -10.0348   -9.8702   -9.9538
   13.7327   13.6474   13.8789   13.9325   13.7827   13.7665   13.7762   13.8681   13.8867   13.7211   14.0961   13.8706
    2.8000    4.3598    0.0020    2.8127    3.5490    4.4005    4.1178    3.1338    0.3334    2.7831    1.8745    1.7856
   15.6977   15.6349   15.7098   15.6858   15.6061   15.6339   15.7103   15.6193   15.6295   15.7294   15.6074   15.5958
    0.0535   -0.2251   -0.1516   -0.0709   -0.2109   -0.1276   -0.0340   -0.1752   -0.0117    0.0436   -0.2481   -0.1603
    1.5652    1.6167    1.4897    1.5206    1.6080    1.6259    1.6008    1.5196    1.7133    1.5597    1.4531    1.7096
    8.3976    8.6141    8.6149    8.5123    8.5710    8.5146    8.4985    8.5494    8.3963    8.4382    8.6104    8.5097
    3.6442    3.2731    3.2538    3.3828    3.2633    3.2936    3.2387    3.4129    3.4920    3.3685    3.5169    3.3055
    8.2565    0.0542    4.5196   -2.2731   -1.6399   -1.7370   -2.9377    0.6546   -1.7598    3.9814    6.7136   -0.3827
    2.5533    2.8328    2.8292    2.7929    2.9282    2.8804    2.9859    2.9519    2.8087    2.7650    2.7130    2.9940
    3.9583    4.1863    4.1261    4.0012    4.1922    3.9954    4.0318    3.7330    3.8616    4.1038    3.8526    3.8081];


ParamColumnToUse = 2;
paramVecToUse = EstimatedParams(:, ParamColumnToUse);
indicesMasterVectorEstimated = [1 3 5 21 23 15 30 31 32 2 6 28 33];

activeNames = {... % param name, nominal value, rage of parameters for uniform prior,
    'TX_elong_glob'                      , exp(2.6),    [exp(0) exp(5)]          %1 from est params above
    'TL_elong_glob'                      , exp(3.5),    [exp(0) exp(6)]          %2 from est params above
    'AGTPdeg_time'                       , exp(8.8),    [exp(6) exp(11)]         %3 from est params above
    'AGTPreg_ON'                         , exp(-3.9),   [exp(-6) exp(-1)]       %4 fixed in mcmc_info_vnprl_F2
    'AGTPdeg_rate'                       , exp(-9.9),   [exp(-13) exp(-7)]       %5 from est params above
    'TXTL_UTR_UTR1_Kd'                   , exp(11),     [exp(-3) exp(15)]         %6 from est params above 
    'TXTL_PTET_RNAPbound_Kd'             , exp(14),     [exp(0) exp(17)]          %7 TO BE ESTIMATED HERE
    'TXTL_PTET_RNAPbound_F'              , exp(1.5),    [exp(0) exp(4)]          %8 fixed in mcmc_info_vnprl_F2
    'TXTL_NTP_RNAP_1_Kd'                 , exp(2.9),    [exp(0) exp(5)]          %9 fixed in mcmc_info_vnprl_F2
    'TXTL_NTP_RNAP_2_Kd'                 , exp(14),     [exp(10) exp(20)]         %10 fixed in mcmc_info_vnprl_F2
    'TXTL_PTET_sequestration_Kd'         , exp(12),     [exp(3) exp(15)]          %11 TO BE ESTIMATED HERE
    'TXTL_PTET_sequestration_F'          , exp(1.5),    [exp(-2) exp(5)]         %12 TO BE ESTIMATED HERE
    'TL_AA_Kd'                           , exp(6.6),    [exp(3) exp(10)]         %13 fixed in mcmc_info_vnprl_F2
    'TL_AGTP_Kd'                         , exp(14.5),   [exp(10) exp(18)]       %14 fixed in mcmc_info_vnprl_F2
    'TXTL_RNAdeg_Kd'                     , exp(15.2),   [exp(7) exp(17)]        %15 from est params above
    'TXTL_INDUCER_TETR_ATC_Kd'           , exp(13),     [exp(0) exp(18)]          %16 TO BE ESTIMATED HERE
    'TXTL_INDUCER_TETR_ATC_F'            , exp(2.6),    [exp(-2) exp(5)]         %17 TO BE ESTIMATED HERE
    'TXTL_DIMER_tetR_Kd'                 , exp(13),     [exp(7) exp(17)]          %18 TO BE ESTIMATED HERE
    'TXTL_DIMER_tetR_F'                  , exp(2.6),    [exp(-2) exp(5)]         %19 TO BE ESTIMATED HERE
    'TXTL_UTR_UTR1_F'                    , exp(-.2),    [exp(-4) exp(2)]         %20 fixed in mcmc_info_vnprl_F2
    'TXTL_PLAC_RNAPbound_Kd'             , exp(13.8),   [exp(5) exp(17)]        %21 from est params above
    'TXTL_PLAC_RNAPbound_F'              , exp(2.6),    [exp(-2) exp(5)]         %22 fixed in mcmc_info_vnprl_F2
    'TXTL_RNAPBOUND_TERMINATION_RATE'    , exp(1.8),    [exp(0) exp(12)]         %23 from est params above
    'TXTL_NTP_RNAP_1_F'                  , exp(0),      [exp(-2) exp(3)]           %24 fixed in mcmc_info_vnprl_F2
    'TXTL_NTP_RNAP_2_F'                  , exp(0),      [exp(-2) exp(3)]           %25 fixed in mcmc_info_vnprl_F2
    'TL_AA_F'                            , exp(-0.3),   [exp(-3) exp(3)]        %26 fixed in mcmc_info_vnprl_F2
    'TL_AGTP_F'                          , exp(-1.2),   [exp(-4) exp(2)]        %27 fixed in mcmc_info_vnprl_F2
    'TXTL_RIBOBOUND_TERMINATION_RATE'    , exp(2.3),    [exp(0) exp(12)]          %28 from est params above
    'TXTL_RNAdeg_F'                      , exp(0),      [exp(-3) exp(3)]           %29 fixed in mcmc_info_vnprl_F2
    'TXTL_RNAdeg_kc'                     , exp(-0.45),   [exp(-5) exp(3)]       %30 from est params above
    'RNAP'                               , exp(1.4419),  [exp(-1) exp(8)]       %31 31% from est params above
    'RNase'                              , exp(8.5),  [exp(5) exp(10)]          %32 from est params above
    'Ribo'                               , exp(3.75),  [exp(1) exp(6)]          %33 from est params above
    'TXTL_PROT_deGFP_MATURATION'         , exp(-6.07), [exp(-9) exp(-3)]    %34 fixed in mcmc_info_vnprl_F2
    'TXTL_INDUCER_LASR_AHL_Kd'           , exp(-2), [exp(5) exp(20)] %36-1
    'TXTL_INDUCER_LASR_AHL_F'            , exp(0), [exp(-6) exp(6)] %37-1
    'TXTL_PLAS_RNAPbound_Kd'             , exp(30), [exp(25) exp(40)]%38-1
    'TXTL_PLAS_RNAPbound_F'              , exp(3), [exp(-6) exp(6)]%39-1 pol_{F,las}
    'TXTL_PLAS_TFBIND_Kd'                , exp(5), [exp(0) exp(10)]%40-1 plas_{tf, Kd} <
    'TXTL_PLAS_TFRNAPbound_Kd'           , exp(8), [exp(0) exp(15)]%41-1 plas-pol_{tf, Kd} <
    'TXTL_PLAS_TFRNAPbound_F'            , exp(0), [exp(-6) exp(6)]%42-1 plas-pol_{tf, F}
    'TXTL_PLAS_TFBIND_F'                 , exp(0), [exp(-6) exp(6)]};%43-1 plas_{tf, F}
% Set the master vector values that are set from the values estimated in "vnprl_F2"
activeNames(indicesMasterVectorEstimated, 2) = num2cell(exp(paramVecToUse));

% Set the master vector values that were already fixed in "vnprl_F2"
preFixedParams = {...
4     'AGTPreg_ON'                         , exp( -3.9120)
34    'TXTL_PROT_deGFP_MATURATION'         , exp( -6.0748)
8     'TXTL_PTET_RNAPbound_F'              , exp( 1.5000)
9     'TXTL_NTP_RNAP_1_Kd'                 , exp( 2.9459)
10    'TXTL_NTP_RNAP_2_Kd'                 , exp( 13.9970)
13    'TL_AA_Kd'                           , exp( 6.5566)
14    'TL_AGTP_Kd'                         , exp( 14.5090)
20    'TXTL_UTR_UTR1_F'                    , exp( -0.2000)
22    'TXTL_PLAC_RNAPbound_F'              , exp( 1.5000)
24    'TXTL_NTP_RNAP_1_F'                  , exp(      0)
25    'TXTL_NTP_RNAP_2_F'                  , exp(      0)
26    'TL_AA_F'                            , exp( -0.3000)
27    'TL_AGTP_F'                          , exp( -1.2000)
29    'TXTL_RNAdeg_F'                      , exp(      0)}; % checked and verified. 

% set the prefixed params elements in master vector to the values in
% prefixed params. 
% 
activeNames(cell2mat(preFixedParams(:,1)),2) = preFixedParams(:,3);

% format short g; 
% aa = [log(cell2mat(activeNames(cell2mat(preFixedParams(:,1)),2))) log(cell2mat(activeNames(cell2mat(preFixedParams(:,1)),3)))]
% 
% pdiagnostic = cell2table([num2cell(aa) preFixedParams(:,1) num2cell([-aa(:,1)+aa(:,3) ...
%     aa(:,1)-aa(:,2) (1:length(aa))']) (activeNames(cell2mat(preFixedParams(:,1)),1))],...
%     'VariableNames', {'logval', 'loglb', 'logub', 'mvix', 'ubdiff', 'lbdiff', 'ix', 'name'})
% pdiagnostic.lbdiff./pdiagnostic.ubdiff
% pdiagnostic =
% 
%   14×8 table
% 
%     logval     loglb    logub    mvix    ubdiff    lbdiff    ix                name            
%     _______    _____    _____    ____    ______    ______    __    ____________________________
% 
%      -3.912     -6       -1        4      2.912     2.088     1    'AGTPreg_ON'                
%     -6.0748     -9       -3       34     3.0748    2.9252     2    'TXTL_PROT_deGFP_MATURATION'
%         1.5      0        4        8        2.5       1.5     3    'TXTL_PTET_RNAPbound_F'     
%      2.9459      0        5        9     2.0541    2.9459     4    'TXTL_NTP_RNAP_1_Kd'        
%      13.997     10       20       10      6.003     3.997     5    'TXTL_NTP_RNAP_2_Kd'        
%      6.5566      3       10       13     3.4434    3.5566     6    'TL_AA_Kd'                  
%      14.509     10       18       14      3.491     4.509     7    'TL_AGTP_Kd'                
%        -0.2     -4        2       20        2.2       3.8     8    'TXTL_UTR_UTR1_F'           
%         1.5     -2        5       22        3.5       3.5     9    'TXTL_PLAC_RNAPbound_F'     
%           0     -2        3       24          3         2    10    'TXTL_NTP_RNAP_1_F'         
%           0     -2        3       25          3         2    11    'TXTL_NTP_RNAP_2_F'         
%        -0.3     -3        3       26        3.3       2.7    12    'TL_AA_F'                   
%        -1.2     -4        2       27        3.2       2.8    13    'TL_AGTP_F'                 
%           0     -3        3       29          3         3    14    'TXTL_RNAdeg_F'    


%% Next set the parameters estimated in mcmc_info_ZSIFFL_mtet_phase1.m
%
% The values estimated are: 
mtet_phase1_params = ...
{...
7		, 	    'TXTL_PTET_RNAPbound_Kd'    	, 		exp(14)             , 		[exp(0)     exp(17)]
11		, 	    'TXTL_PTET_sequestration_Kd'	, 		exp(-1)             , 		[exp(-10)	exp(5)]
12		, 	    'TXTL_PTET_sequestration_F' 	, 		exp(1.314)          , 		[exp(-2)	exp(5)]
16		, 	    'TXTL_INDUCER_TETR_ATC_Kd'  	, 		exp(-2)             , 		[exp(-15)	exp(5)]
17		, 	    'TXTL_INDUCER_TETR_ATC_F'   	, 		exp(1.577)          , 		[exp(-2)	exp(5)]
18		, 	    'TXTL_DIMER_tetR_Kd'        	, 		exp(-10)            , 		[exp(-20)	exp(-7)]
19		, 	    'TXTL_DIMER_tetR_F'         	, 		exp(1.447)          , 		[exp(-2)	exp(5)]...
};
activeNames(cell2mat(mtet_phase1_params(:,1)),2) = mtet_phase1_params(:,3);
activeNames(cell2mat(mtet_phase1_params(:,1)),3) = mtet_phase1_params(:,4);

%% next we set the parameters estimated from mcmc_info_ZSIFFL_mtet_phase2.m

% The values estimated are: 
mtet_phase2_params = ...
{...
1       ,       'TX_elong_glob'                 , exp(2.3)      ,[exp(0)    exp(5)]          %1 from est params above
2       ,       'TL_elong_glob'                 , exp(3.7)      ,[exp(0)    exp(6)]          %2 from est params above
3       ,       'AGTPdeg_time'                  , exp(10.05)    ,[exp(8)    exp(12)]         %3 from est params above
31      ,       'RNAP'                          , exp(5.9)      ,[exp(1)   exp(15)]       %31 31% from est params above
32      ,       'RNase'                         , exp(9.2)      ,[exp(7)    exp(11)]          %32 from est params above
33      ,       'Ribo'                          , exp(5.9)      ,[exp(1)    exp(15)]          %33 from est params above
21      ,       'TXTL_PLAC_RNAPbound_Kd'        , exp(10.5)     ,[exp(5)    exp(17)]        %21 from est params above
7		, 	    'TXTL_PTET_RNAPbound_Kd'    	, exp(17)    ,[exp(0)    exp(22)]
11		, 	    'TXTL_PTET_sequestration_Kd'	, exp(-2.7)     ,[exp(-10)	exp(5)]
16		, 	    'TXTL_INDUCER_TETR_ATC_Kd'  	, exp(-6)       ,[exp(-15)	exp(5)]...
};
activeNames(cell2mat(mtet_phase2_params(:,1)),2) = mtet_phase2_params(:,3);
activeNames(cell2mat(mtet_phase2_params(:,1)),3) = mtet_phase2_params(:,4);

%% we also set the forward rate parameters, since those should not matter much, and any value around 1 is good
% the F rate parameters just set the timescale. 

mtet_phase2_params = ...
{...
36      ,       'TXTL_INDUCER_LASR_AHL_F'     	,   exp(0)        , [exp(-3)    exp(5)] 
38      ,       'TXTL_PLAS_RNAPbound_F'         ,   exp(0)        , [exp(-3)    exp(5)] 
41      ,       'TXTL_PLAS_TFRNAPbound_F'       ,   exp(0)        , [exp(-3)    exp(5)] 
42      ,       'TXTL_PLAS_TFBIND_F'            ,   exp(0)        , [exp(-3)    exp(5)] };
activeNames(cell2mat(mtet_phase2_params(:,1)),2) = mtet_phase2_params(:,3);
activeNames(cell2mat(mtet_phase2_params(:,1)),3) = mtet_phase2_params(:,4);



estParamsIX = [23 28 31 33 35 37 39 40]';
estParams = activeNames(estParamsIX,1);
activeNames(estParamsIX,:);
% skipping AGTPdeg_rate, AGTPreg_ON, TXTL_PROT_deGFP_MATURATION
% fixedParams vector
fixedParamsIX =  setdiff((1:size(activeNames, 1))', estParamsIX);

% for troubleshooting / visualizing the fixed params:
%     log(cell2mat(activeNames2(fixedParamsIX,[2])))
%     activeNames2(fixedParamsIX,[1:2])


% since activeNames2 is a superset of activeNames1, we can just use
% activeNames2 as the master vector.
masterVector = log(cell2mat(activeNames(:,2))); % log transformed.

% paramMap is a matrix mapping the parameters in the master vector to the
% (unordered) list of parameters in the model. (obvioulsy within the code
% these parameters get ordered before they are used in the exported model)
% More precisely, Let pp = paramMap(:, 1); then masterVector(pp) is the list
% of parameters for the first geometry within that topology, as specified by
% namesUnord. Note that namesUnord is just all the active parameters in
% the model, not just the estimated ones.
% One such matrix exists for each topology. It has dimnesions
% length(model_info(i).namesUnord) x number of geometries associated with that topo.

paramMap_ptet = [1:10 13:15 20 23:34]';
paramMap_tetR = [1:15 18:34]';
paramMap_aTc = (1:34)';
paramMap_plac = [1:6 9 10 13:15 20:34]';
paramMap_plas = [1:6 9 10 13:15 20:42]';

%     'TX_elong_glob'                      , exp(2.6),   [exp(0) exp(5)]          %1 from est params above
%     'TL_elong_glob'                      , exp(3.5),   [exp(0) exp(6)]          %2 from est params above
%     'AGTPdeg_time'                       , exp(8.8),   [exp(6) exp(11)]         %3 from est params above
%     'AGTPreg_ON'                         , exp(-3.9),   [exp(-6) exp(-1)]       %4 fixed in mcmc_info_vnprl_F2
%     'AGTPdeg_rate'                       , exp(-9.9),  [exp(-13) exp(-7)]       %5 from est params above
%     'TXTL_UTR_UTR1_Kd'                   , exp(11),   [exp(-3) exp(15)]         %6 from est params above 
%     'TXTL_NTP_RNAP_1_Kd'                 , exp(2.9),   [exp(0) exp(5)]          %9 fixed in mcmc_info_vnprl_F2
%     'TXTL_NTP_RNAP_2_Kd'                 , exp(14),   [exp(10) exp(20)]         %10 fixed in mcmc_info_vnprl_F2
%     'TL_AA_Kd'                           , exp(6.6),   [exp(3) exp(10)]         %13 fixed in mcmc_info_vnprl_F2
%     'TL_AGTP_Kd'                         , exp(14.5),   [exp(10) exp(18)]       %14 fixed in mcmc_info_vnprl_F2
%     'TXTL_RNAdeg_Kd'                     , exp(15.2),   [exp(7) exp(17)]        %15 from est params above
%     'TXTL_UTR_UTR1_F'                    , exp(-.2),   [exp(-4) exp(2)]         %20 fixed in mcmc_info_vnprl_F2
%     'TXTL_PLAC_RNAPbound_Kd'             , exp(13.8),   [exp(5) exp(17)]        %21 from est params above
%     'TXTL_PLAC_RNAPbound_F'              , exp(2.6),   [exp(-2) exp(5)]         %22 fixed in mcmc_info_vnprl_F2
%     'TXTL_RNAPBOUND_TERMINATION_RATE'    , exp(1.8),   [exp(-3) exp(6)]         %23 from est params above
%     'TXTL_NTP_RNAP_1_F'                  , exp(0),   [exp(-2) exp(3)]           %24 fixed in mcmc_info_vnprl_F2
%     'TXTL_NTP_RNAP_2_F'                  , exp(0),   [exp(-2) exp(3)]           %25 fixed in mcmc_info_vnprl_F2
%     'TL_AA_F'                            , exp(-0.3),   [exp(-3) exp(3)]        %26 fixed in mcmc_info_vnprl_F2
%     'TL_AGTP_F'                          , exp(-1.2),   [exp(-4) exp(2)]        %27 fixed in mcmc_info_vnprl_F2
%     'TXTL_RIBOBOUND_TERMINATION_RATE'    , exp(2.3),   [exp(0) exp(5)]          %28 from est params above
%     'TXTL_RNAdeg_F'                      , exp(0),   [exp(-3) exp(3)]           %29 fixed in mcmc_info_vnprl_F2
%     'TXTL_RNAdeg_kc'                     , exp(-0.45),   [exp(-5) exp(3)]       %30 from est params above
%     'RNAP'                               , exp(1.4419),  [exp(-1) exp(4)]       %31 31% from est params above
%     'RNase'                              , exp(8.5),  [exp(5) exp(10)]          %32 from est params above
%     'Ribo'                               , exp(3.75),  [exp(1) exp(6)]          %33 from est params above
%     'TXTL_PROT_deGFP_MATURATION'         , exp(-6.07), [exp(-9) exp(-3)]  };    %34 fixed in mcmc_info_vnprl_F2

    
% parameter ranges (for the to-be-estimated parameters in the master
% vector)
paramRanges = log(cell2mat(activeNames(estParamsIX,3)));


%% next we define the dosing strategy.

dosedNames1 = {'DNA ptet--utr1--deGFP'};
dosedVals1 = [4 2 1 0.5 0.25 0.125 0.0625];
%     dtempvec = sqrt(dosedVals1(end)*(ones(size(dosedVals1))./dosedVals1)).*dosedVals1;
doseWeights1 = ones(1,size(dosedVals1,2)); %dtempvec/(sum(dtempvec));

dosedNames2 = {'DNA ptet--utr1--deGFP'; 
    'DNA plac--utr1--tetR'};
dosedVals2 = [1 1 1 1 1 1 1; 
    2 0.2 0.02 0.002 0.0002 0.00002 0.000002];
%     dtempvec = sqrt(dosedVals2(end)*(ones(size(dosedVals2))./dosedVals2)).*dosedVals2;
doseWeights2 = ones(1,size(dosedVals2,2)); %dtempvec/(sum(dtempvec));

dosedNames3 = {'DNA ptet--utr1--deGFP'; 
    'DNA plac--utr1--tetR';
    'aTc'};
dosedVals3 = [1 1 1 1 1 1 1;
    0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1;
    10000 1000 100 10 1 0.1 0.01];
doseWeights3 = ones(1,size(dosedVals3,2)); %dtempvec/(sum(dtempvec));

dosedNames4 = {'DNA plac--utr1--deGFP'};
dosedVals4 = [2 1 0.5 0.25 0.125 0.0625 0.0313];

doseWeights4 = ones(1,size(dosedVals4,2)); %dtempvec/(sum(dtempvec));


dosedNames5 = {...
    'DNA plac--utr1--lasR';
    'DNA plas--utr1--deGFP';
    'OC12HSL'};
dosedVals5 = [1 1 1 1 1 1 1;
    1 1 1 1 1 1 1;
    10000 , 1000, 100, 10, 1, 0.1, 0.01];

doseWeights5 = ones(1,size(dosedVals5,2)); %dtempvec/(sum(dtempvec));




%% create the measured species cell array
% remember to change this! esp the 2AGTP.
measuredSpecies = {{'protein deGFP*'}};
msIx = 1; % this is the index of the measured species in the data array


%% setup the MCMC simulation parameters
stdev = 100; % i have no idea what a good value is
tightening = 1; % i have no idea what a good value is
nW = 30; % actual: 200 - 600 ish
stepsize = 1.2; % actual: 1.1 to 4 ish
niter =3; % actual: 2 - 30 ish,
npoints = 1e3; % actual: 2e4 to 2e5 ish (or even 1e6 of the number of
%                        params is small)
thinning = 10; % actual: 10 to 40 ish

%% pull all this together into an output struct.

runsim_info = struct('stdev', {stdev}, ...
    'tightening', {tightening}, ...
    'nW', {nW}, ...
    'stepSize', {stepsize}, ...
    'nIter', {niter}, ...
    'nPoints', {npoints}, ...
    'thinning', {thinning}, ...
    'parallel', false);

mtet1 = copyobj(mtet);
mtet2 = copyobj(mtet);
mtet3 = copyobj(mtet);
model_info = struct(...
    'circuitInfo',{ptetinfo, tetrinfo, atcinfo, placinfo, plasinfo},...
    'modelObj', {mtet1, mtet2, mtet3, mlac, mlas},... % array of model objects (different topologies)
    'modelName',  {mtet1.name, mtet2.name, mtet3.name, mlac.name, mlas.name},...; % model names.
    'namesUnord', {activeNames(paramMap_ptet,1),activeNames(paramMap_tetR,1),activeNames(paramMap_aTc,1),activeNames(paramMap_plac,1), activeNames(paramMap_plas,1)},... % names of parameters per model, unordered.
    'paramMaps', {paramMap_ptet, paramMap_tetR, paramMap_aTc, paramMap_plac, paramMap_plas}, ... % paramMap is a matrix mapping master vector elements to namesUnord
    'dosedNames', {dosedNames1, dosedNames2, dosedNames3, dosedNames4, dosedNames5},... % cell arrays of species. cell array corresponds
    ...                               % to a model.
    'dosedVals', {dosedVals1, dosedVals2, dosedVals3, dosedVals4, dosedVals5},...  % matrices of dose vals
    'doseWeighting', {doseWeights1, doseWeights2, doseWeights3, doseWeights4, doseWeights5}, ...
    ... % OPTIONAL FIELD. reweight the importance of the curves corresponding to the different doses.
    'measuredSpecies', {measuredSpecies, measuredSpecies, measuredSpecies , measuredSpecies, measuredSpecies}, ... % cell array of cell arrays of
    ...                  % species names. the elements of the inner
    ...                  % cell array get summed.
    'measuredSpeciesIndex', {msIx, msIx,msIx, msIx, msIx},...  % maps measuredSpecies to the species in data array
    'experimentWeighting', {1, 1, 1, 1, 1}, ... %
    ... % relative importance of the different topologies.
    ... %geometries in a given topology are weighted with
    ...% the same level of importance for now.
    'dataToMapTo', {1,3,4, 2, 5}); % each dataToMapTo property within an element of the
% model_info array is a vector of length # of geometries.
% data indices tell us which data set to use for each topology (model) - geometry pair
% from the data_info struct array.

semanticGroups = num2cell((1:length(estParams))');
%arrayfun(@num2str, 1:10, 'UniformOutput', false);
% estParamsIx = setdiff((1:length(masterVector))', fixedParamsIX);

%% master parameter vector, param ranges,
master_info = struct(...
    'estNames', {estParams},...
    'masterVector', {masterVector},...
    'paramRanges', {paramRanges},... %
    'fixedParams', {fixedParamsIX},...   % indexes of the fixed params (withing master vector)
    'semanticGroups', {semanticGroups}); % EITHER EMPTY OR
% a cell array of vectors specifying parameter
% groupings.
% The vectors contain indices to the
% parameters in (non fixed subset of) the master
% vector that need to be grouped.
% I.e., They contain indexes of the subvector
% logp =
% master_info.mastervector(~master_info.fixedParams)
% and to the rows of the paramRanges matrix and the
% estNames cell array of strings.
%
% parameter grouping so that these parameters
% get INITIALIZED to the same values.
%
% every parameter index must show up in at least
% one group, even if that is the only parameter in
% that group. If the semanticGroups field is empty,
% then all parameters are assumed to be in their
% distinct groups.


% how the parameter distribution flow works:
% WALKER INITIALIZATION
% reduced master vector -- semanticGroups -->
% master vector -- paramMaps -->
% full parameter vector for each topo-geom pair -- orderingIx -->
% reordered vector for exported model simulation.
%
% param ranges: reduced param ranges matrix (by sematicGroups)
% compute initial parameter distributons
% then expand in the same way as above.
% once the parameters have been estimated, there is no need to
% reorder them, since the master vector was never reordered.
% can use the master_info.estNames for the names and
% master_info.mastervector(~master_info.fixedParams) for the
% parameter values.
% sematic


mcmc_info = struct('runsim_info', runsim_info, ...
    'model_info', model_info,...
    'master_info', master_info);


if nargout == 2
    varargout{1} = activeNames;
    
end

end