function [mcmc_info, varargout] = mcmc_info_ZSIFFL_mtet(mtet)
% This project sets up the concurrent estimation problem for the IFFL data
% generated by Zach Sun in 2014. 
%
% plac - UTR1 - tetR, ptet - UTR1 - deGFP, aTc; This has three
% topologies, the constitutive expression geometry, the repression
% geometry, and the aTc induction geometry. 

% 
% The parameters are shared between these five circuits (experiments). 
%
% A lot of the core parameters are set from those estimated from the
% constitutive expression fits to the data in the VNPRL2011 and ACS 2014 papers. 
%
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
    
    % first we set up the master vector
    % 
    % then we will set the parameters to fix, and the remaining ones to
    % estimate. 
    %
    % Then we will define paramMaps that will be used to distribute the
    % master vector among the different topologies. 
    
    
    
    masterVector = {'TX_elong_glob'  , '10.5'          
'TL_elong_glob'                      , '20'            
'AGTPdeg_time'                       , '7200'          
'AGTPreg_ON'                         , '0.02'          
'AGTPdeg_rate'                       , '0.0002'
'TXTL_PROT_deGFP_MATURATION'         , '0.00231049' 
'TXTL_UTR_UTR1_Kd'                   , '20'            
'TXTL_PTET_RNAPbound_Kd'             , '20'            
'TXTL_PTET_RNAPbound_F'              , '20'            
'TXTL_NTP_RNAP_1_Kd'                 , '100000'        
'TXTL_NTP_RNAP_2_Kd'                 , '1e+06'         
'TXTL_PTET_sequestration_Kd'         , '0.000444'      
'TXTL_PTET_sequestration_F'          , '0.25'          
'TL_AA_Kd'                           , '100000'        
'TL_AGTP_Kd'                         , '100000'        
'TXTL_RNAdeg_Kd'                     , '200'           
'TXTL_INDUCER_TETR_ATC_Kd'           , '10000'         
'TXTL_INDUCER_TETR_ATC_F'            , '1e-05'          
'TXTL_DIMER_tetR_Kd'                 , '1'             
'TXTL_DIMER_tetR_F'                  , '1'              
'TXTL_UTR_UTR1_F'                    , '0.2'          
'TXTL_PLAC_RNAPbound_Kd'             , '20'            
'TXTL_PLAC_RNAPbound_F'              , '20'            
'TXTL_RNAPBOUND_TERMINATION_RATE'    , '0.15'          
'TXTL_NTP_RNAP_1_F'                  , '0.0001'         
'TXTL_NTP_RNAP_2_F'                  , '1e-05'           
'TL_AA_F'                            , '0.001'          
'TL_AGTP_F'                          , '1e-05'          
'TXTL_RIBOBOUND_TERMINATION_RATE'    , '40'            
'TXTL_RNAdeg_F'                      , '0.01'           
'TXTL_RNAdeg_kc'                     , '0.00277778'  };
    
    
    % parameters of the master vector we fix vs estimate
    
masterVector = {'TX_elong_glob'      , '10.5'          
'TL_elong_glob'                      , '20'            
'AGTPdeg_time'                       , '7200'          
'AGTPreg_ON'                         , '0.02'          
'AGTPdeg_rate'                       , '0.0002'        
'TXTL_UTR_UTR1_Kd'                   , '20'            
'TXTL_PTET_RNAPbound_Kd'             , '20'            
'TXTL_PTET_RNAPbound_F'              , '20'            
'TXTL_NTP_RNAP_1_Kd'                 , '100000'        
'TXTL_NTP_RNAP_2_Kd'                 , '1e+06'         
'TXTL_PTET_sequestration_Kd'         , '0.000444'      
'TXTL_PTET_sequestration_F'          , '0.25'          
'TL_AA_Kd'                           , '100000'        
'TL_AGTP_Kd'                         , '100000'        
'TXTL_RNAdeg_Kd'                     , '200'           
'TXTL_INDUCER_TETR_ATC_Kd'           , '10000'         
'TXTL_INDUCER_TETR_ATC_F'            , '1e-05'          
'TXTL_DIMER_tetR_Kd'                 , '1'             
'TXTL_DIMER_tetR_F'                  , '1'              
'TXTL_UTR_UTR1_F'                    , '0.2'          
'TXTL_PLAC_RNAPbound_Kd'             , '20'            
'TXTL_PLAC_RNAPbound_F'              , '20'            
'TXTL_RNAPBOUND_TERMINATION_RATE'    , '0.15'          
'TXTL_NTP_RNAP_1_F'                  , '0.0001'         
'TXTL_NTP_RNAP_2_F'                  , '1e-05'           
'TL_AA_F'                            , '0.001'          
'TL_AGTP_F'                          , '1e-05'          
'TXTL_RIBOBOUND_TERMINATION_RATE'    , '40'            
'TXTL_RNAdeg_F'                      , '0.01'           
'TXTL_RNAdeg_kc'                     , '0.00277778'  };
    
    
    
    
    activeNames1 = {...
        'TXTL_RNAdeg_Kd'                    exp(9.237)        [100 1e5] 
        'TXTL_RNAdeg_F'                     exp(0)        [0.01 10000]
        'TXTL_RNAdeg_kc'                    exp(-4.4)      [1e-4 1]
        'RNase'                             exp(6.4899)         [1 10000]};
    activeNames2 = {... % changes made to ranges on feb 8, 2019. setting parameters based on 
        ...% posterior plots. 
        'TX_elong_glob'                     exp(4.9354)       [0.5 300]        % 1
        'AGTPdeg_time'                      exp(9.17)   [1800 42000] % DO NOT FIX THIS. 
        'AGTPdeg_rate'                      exp(-9.5172)      [1e-7 1e-2] % set from before
        'AGTPreg_ON'                        exp(-3.912)        [0.005 0.2]        %4 % set from before
        'TXTL_P70_RNAPbound_Kd'             exp(9.514)         [0.01 1e7]  % 
        'TXTL_P70_RNAPbound_F'              exp(1.5)    [1e-5 300] % set to exp(1.5)
        'TXTL_RNAPBOUND_TERMINATION_RATE'   exp(3.3005)        [0.1 1000]         % 7
        'TXTL_NTP_RNAP_1_Kd'                exp(2.9459)      [1 1e6]
        'TXTL_NTP_RNAP_1_F'                 exp(0)      [1e-5 100] % set to 1
        'TXTL_NTP_RNAP_2_Kd'                exp(13.997)         [0.1 1e7]
        'TXTL_NTP_RNAP_2_F'                 exp(0)      [1e-6 1000]        %11 % set to 1
        'TXTL_RNAdeg_Kd'                    exp(9.237)        [100 1e7] 
        'TXTL_RNAdeg_F'                     exp(0)        [0.01 10000] % set to 1 (1 is right in the middle of the broad posterior density, and so not entirely arbitrary. )
        'TXTL_RNAdeg_kc'                    exp(-4.4)   [1e-5 10]  %set to exp(-5.4)
        'RNAP'                              exp(1.4419)         [0.1 5000] % 15
        'RNase'                             exp(6.4899)         [1 10000]
        'TL_elong_glob'                     exp(0.5219)          [0.1 500]
        'TXTL_PROT_deGFP_MATURATION'        exp(-6.0748)      [0.0002 0.02] %18 % set from before
        'TXTL_UTR_UTR1_Kd'                  exp(11.189)          [0.05 1e7]
        'TXTL_UTR_UTR1_F'                   exp(-0.2)   [1e-5 100] % set to exp(-0.2)
        'TL_AA_Kd'                          exp(6.5566)      [.1 1e6] % 21
        'TL_AA_F'                           exp(-0.3)   [1e-5 20]   % set to exp(-0.3)
        'TL_AGTP_Kd'                        exp(14.509)      [.1 1e7] % 23
        'TL_AGTP_F'                         exp(-1.2)   [1e-5 100]  %set to exp(-1.2)
        'TXTL_RIBOBOUND_TERMINATION_RATE'   exp(5.398)          [0.1 1000]
        'Ribo'                              exp(7.3081)          [10 10000]}; %26 
    
    %%
    % Names of parameters and species to actually estimate.
    estParamsIX = [1 2 3 5 7 12 14 15 16 17 19 25 26]';
    estParams = activeNames2(estParamsIX,1);
    % skipping AGTPdeg_rate, AGTPreg_ON, TXTL_PROT_deGFP_MATURATION
    % fixedParams vector
    fixedParamsIX =  setdiff((1:26)', estParamsIX);
    
    % for troubleshooting / visualizing the fixed params:
%     log(cell2mat(activeNames2(fixedParamsIX,[2])))
%     activeNames2(fixedParamsIX,[1:2])


    % since activeNames2 is a superset of activeNames1, we can just use
    % activeNames2 as the master vector.
    masterVector = log(cell2mat(activeNames2(:,2))); % log transformed.
    
    % paramMap is a matrix mapping the parameters in the master vector to the
    % (unordered) list of parameters in the model. (obvioulsy within the code
    % these parameters get ordered before they are used in the exported model)
    % More precisely, Let pp = paramMap(:, 1); then masterVector(pp) is the list
    % of parameters for the first geometry within that topology, as specified by
    % namesUnord. Note that namesUnord is just all the active parameters in
    % the model, not just the estimated ones.
    % One such matrix exists for each topology. It has dimnesions
    % length(model_info(i).namesUnord) x number of geometries associated with that topo.
    paramMap1 = [12 13 14 16]';
    paramMap2 = (1:length(masterVector))';
    
    
    % parameter ranges (for the to-be-estimated parameters in the master
    % vector)
    paramRanges = log(cell2mat(activeNames2(estParamsIX,3)));
    
    
    %% next we define the dosing strategy.
    
    dosedNames1 = {'RNA utr1--deGFP'};
    dosedVals1 = [37.5 75 150 200 600 700 800 900 1000];
%     dtempvec = sqrt(dosedVals1(end)*(ones(size(dosedVals1))./dosedVals1)).*dosedVals1;
    doseWeights1 = ones(size(dosedVals1)); %dtempvec/(sum(dtempvec));
    dosedNames2 = {'DNA p70--utr1--deGFP'};
    dosedVals2 = [0.5 2 5 20];
%     dtempvec = sqrt(dosedVals2(end)*(ones(size(dosedVals2))./dosedVals2)).*dosedVals2;
    doseWeights2 = ones(size(dosedVals2)); %dtempvec/(sum(dtempvec));
    %% create the measured species cell array
    % remember to change this! esp the 2AGTP.
    measuredSpecies1 = {{'[RNA utr1--deGFP]',...
        '[Ribo:RNA utr1--deGFP]',...
        '[AA:AGTP:Ribo:RNA utr1--deGFP]', ...
        '[term_Ribo:RNA utr1--deGFP]',...
        '[AA:Ribo:RNA utr1--deGFP]'...
        '[RNA utr1--deGFP:RNase]',...
        '[Ribo:RNA utr1--deGFP:RNase]',...
        '[AA:AGTP:Ribo:RNA utr1--deGFP:RNase]', ...
        '[term_Ribo:RNA utr1--deGFP:RNase]',...
        '[AA:Ribo:RNA utr1--deGFP:RNase]'}};
    
    measuredSpecies2 = {{'[RNA utr1--deGFP]',...
        '[Ribo:RNA utr1--deGFP]',...
        '[AA:AGTP:Ribo:RNA utr1--deGFP]', ...
        '[term_Ribo:RNA utr1--deGFP]',...
        '[AA:Ribo:RNA utr1--deGFP]'...
        '[RNA utr1--deGFP:RNase]',...
        '[Ribo:RNA utr1--deGFP:RNase]',...
        '[AA:AGTP:Ribo:RNA utr1--deGFP:RNase]', ...
        '[term_Ribo:RNA utr1--deGFP:RNase]',...
        '[AA:Ribo:RNA utr1--deGFP:RNase]'}, {'protein deGFP*'}};
    msIx1 = 1; % this is the index of the measured species in the data array
    % from data_dsg2014. There are two species: 1: mRNA and 2: GFP.
    msIx2 = [1,2];
    
    %% setup the MCMC simulation parameters
    stdev = 100; % i have no idea what a good value is
    tightening = 1; % i have no idea what a good value is
    nW = 300; % actual: 200 - 600 ish
    stepsize = 1.5; % actual: 1.1 to 4 ish
    niter = 10; % actual: 2 - 30 ish,
    npoints = 1e4; % actual: 2e4 to 2e5 ish (or even 1e6 of the number of
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
    
    model_info = struct(...
        'circuitInfo',{circuitInfo1, circuitInfo2},...
        'modelObj', {modelObj,modelObj},... % array of model objects (different topologies)
        'modelName',  modelObj.name,...; % model names.
        'namesUnord', {activeNames1(:,1),activeNames2(:,1)}, ... % names of parameters per model, unordered.
        'paramMaps', {paramMap1, paramMap2}, ... % paramMap is a matrix mapping master vector elements to namesUnord
        'dosedNames', {dosedNames1, dosedNames2},... % cell arrays of species. cell array corresponds
        ...                               % to a model.
        'dosedVals', {dosedVals1, dosedVals2},...  % matrices of dose vals
        'doseWeighting', {doseWeights1, doseWeights2}, ...
        ... % OPTIONAL FIELD. reweight the importance of the curves corresponding to the different doses.
        'measuredSpecies', {measuredSpecies1, measuredSpecies2}, ... % cell array of cell arrays of
        ...                  % species names. the elements of the inner
        ...                  % cell array get summed.
        'measuredSpeciesIndex', {msIx1, msIx2},...  % maps measuredSpecies to the species in data array
        'experimentWeighting', {1 0.3}, ... % 
        ... % relative importance of the different topologies.
        ... %geometries in a given topology are weighted with
        ...% the same level of importance for now. 
        'dataToMapTo', {3,1}); % each dataToMapTo property within an element of the
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
    
    
    if nargout == 1
        varargout{1} = activeNames2;
        
    end
    
end