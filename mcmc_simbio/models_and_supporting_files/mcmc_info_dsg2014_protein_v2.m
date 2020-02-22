function mcmc_info = mcmc_info_dsg2014_protein_v2(modelObj)
    % version2 mcmc_info struct. Compatible with multimodel parameter inference. 
    % key differences with version 1: 
    % model parameters specified now include the ones to be estimated and the 
    % ones to keep fixed. 
    % split the old mcmc_info struct into 3 structs: 
    % runsim_info:  information on the mcmc algorithm parameters
    % model_info:   array of models, and associated properties like parameters, 
    %               and the matrices of indices from the master vector 
    %               to the model parameters. 
    % master_info:  contains the master vector, and a spec for which parameters
    %               get estimated. 
    % 
    % Version 1's help text: (!TODO: write this version's full help file.)
    % 
    % Define the mcmc_info struct for the dsg2014 dataset, for the estimation 
    % of the protein parameters, using only measurements of the protein for 
    % the estimation. The data is from figure 1 of the paper:
    % Gene Circuit Performance Characterization and Resource Usage in a 
    % Cell-Free “Breadboard” by Siegal-Gaskins et al. 
    %
    % INPUTS: A Simbiology model object. 
    % 
    % OUTPUTS: You should set up an mcmc_info struct with the fields:
    %
    % The strust has fields: 
    % 
    % 'circuitInfo': A human readable despription of the circuit. (optional)
    % 
    % 'modelObj': The simbiology model object
    %
    % 'modelName': The name property of the Simbiology model object. (optional)
    %
    % 'namesUnord': List of species and parameters that are to be estimated. These 
    % are strings naming things in the model object. 
    %
    % 'paramRanges': A length(mcmc_info.namesUnord) x 2 matrix of log transformed 
    % upper and lower bounds for the parameters and species concentrations. 
    %
    % 'dosedNames': A cell array of strings of species names for species that are 
    % dosed in the model. 
    %
    % 'dosedVals': A matrix of dose values of size # of dosed species by 
    % # of dose combinations. 
    %
    % 'measuredSpecies': A 1 x number of measured output variables. This is a 
    % cell array of cell arrays of the strings of species whose concentrations
    %  are to be added to get the measured variable. 
    %
    % 'stdev': MCMC likelihood function standard deviation
    %
    % 'tightening': A division factor for the standard deviation. 
    %
    % 'nW': Number of Walkers
    %
    % 'stepSize': MCMC step size
    %
    % 'nIter': Number of MCMC iterations. 
    %
    % 'nPoints': Number of MCMC points per iteration. 
    %   
    % 'thinning': Number of steps to skip before taking an MCMC sample. 
    %
    % 'parallel': Boolean variable specifying whether parallel computing is used. 
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

    % User readable description of the circuit. Will be used in the log file generated
    % from the MCMC inference procedure. 
    circuitInfo = ...
    ['This is a simple constitutive gene expression model \n'...
    'built using the TXTL modeling toolbox. It models DNA binding \n'...
    'to RNAP and nucleotides, followed by transcription. The resulting\n'...
    'mRNA can degrade and participate in translation. The former is \n'...
    'modeled as a enzymatic reaction involving every complex containing \n'...
    'mRNA. The latter involves binding to Ribosomes, followed by amino acids \n'...
    'and ATP, and finally elongation and termination resulting in protein.']

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
    
    
    
    activeNames = {'TX_elong_glob'
                'AGTPdeg_time'
                'AGTPdeg_rate'
                'TXTL_P70_RNAPbound_Kd'
                'TXTL_P70_RNAPbound_F'
                'TXTL_RNAPBOUND_TERMINATION_RATE'
                'TXTL_NTP_RNAP_1_Kd'
                'TXTL_NTP_RNAP_1_F'
                'TXTL_NTP_RNAP_2_Kd'
                'TXTL_NTP_RNAP_2_F'
                'TXTL_RNAdeg_Kd'
                'TXTL_RNAdeg_F'
                'TXTL_RNAdeg_kc'
                'RNAP'
                'RNase'
                'TL_elong_glob'
                'TXTL_PROT_deGFP_MATURATION'
                'TXTL_UTR_UTR1_Kd'
                'TXTL_UTR_UTR1_F'
                'TL_AA_Kd'
                'TL_AA_F'
                'TL_AGTP_Kd'
                'TL_AGTP_F'
                'TXTL_RIBOBOUND_TERMINATION_RATE'
                'Ribo'};
            
            
estParams = {'TL_elong_glob'
                'TXTL_PROT_deGFP_MATURATION'
                'TXTL_UTR_UTR1_Kd'
                'TXTL_UTR_UTR1_F'
                'TL_AA_Kd'
                'TL_AA_F'
                'TL_AGTP_Kd'
                'TL_AGTP_F'
                'TXTL_RIBOBOUND_TERMINATION_RATE'
                'Ribo'};
            


    % get 10 sets of parameter values (= 10 points) from a previously estimated sims:
    % get parameter values from a few different mcmc points from the runs of 
    % the file proj_acs_dsg2014_protein
    marray = mcmc_get_walkers({'20180121_131114'}, {5}, ...
      ['/Users/vipulsinghal/Dropbox/Documents/toolbox/txtlsim_vsfork2017/'...
      'mcmc_simbio/projects/proj_acs_dsg2014_mrna']);
  
    
    nPrevPoints = 5;
    minit = marray(:, ((end-nPrevPoints+1) : end), end); % Final set of walker positions. 
    minit = minit(:,:); % nparam x npoints, nparam is the tx params, ie, 15. npoints is 10. 

    % fixedParams vector
    fixedParams = [1:numel(minit)]';

    % master vector
    logp =  zeros(10,1);

    masterVector = [minit(:)
                    logp]; % log transformed. 
                    
    % paramMap is a matrix mapping the parameters in the master vector to the 
    % (unordered) list of parameters in the model. (obvioulsy within the code 
    % these parameters get ordered before they are used in the exported model)
    % More precisely, Let pp = paramMap(:, 1); then masterVector(pp) is the list 
    % of parameters for the first geometry within that topology. 
    % 
    % One such matrix exists for each topology. It has dimnesions 
    % length(model_info(i).namesUnord) x number of geometries associated with that topo. 
    % 
    % 

    TLparamIX = [(length(masterVector)-length(logp)+1) : length(masterVector)]';
    paramMap = [reshape((1:numel(minit))', 15, numel(minit)/15); 
                repmat(TLparamIX, 1, numel(minit)/15)];

    % parameter ranges for the logp (the parameters within the master vector that 
    % are to be estimated.)
    
    % 20180218_225205 -> 3 successful runs
    % 20180219_021944 -> 10 successpul runs, comtinued from 20180218_225205
    % used: 
%     paramRanges = [4 8
%     -7 -3
%     -1 7
%     -3 3
%     -1 6
%     -1 6
%     -1 6
%     -1 6
%     -6 -1
%     2 7];

% '20180219_154541', 10 successful runs, started from 
% 20180219_021944, but with the parameter intiial points constrained to:

% paramRanges = [5 8
%     -4 0
%     0 6
%     -1 3
%     0 6
%     0 6
%     -1 6
%     -1 6
%     -1.5 2
%     6 8];


% New run '20180221_011619' started from 20180219_154541 but with the 
% parameter intiial points constrained to:
paramRanges = [3 8
    -4 0
    -2 8
    -2 4
    -2 8
    -2 8
    4 8
    -3 2
    -1.5 -0.5
    3 9];

% New run _______ started from 20180221_011619 but with the 
% parameter intiial points constrained to:
paramRanges = [5 9
    -4 0
    -2 10
    -1 3
    -2 10
    -3 10
    2 7
    -3 4
    -1.2 -0.8
    5 8];


% data indices tell us which data set to use for each topology (model) - geometry pair
% from the data_info struct array. 
% Here we have 1 model, and 10 geometries, all pointing to the same data: 
% at data index 1 in the data_info array created by data_dsg2014.m
dataIndices = ones(nPrevPoints,1);

%% next we define the dosing strategy. 
dosedNames = {'DNA p70--utr1--deGFP'};
dosedVals = [0.5 2 5 20];


%% create the measured species cell array
% this is a 1x2 cell array. each element of this cell array contains
% further cell arrays. The first such cell array is a list of all the bound
% and free versions of the RNA. The second cell array contains a single
% cell, which contains the GFP string. 
% all the species in the inner cell arrays get summed, and compared to the
% corresponding column (dimension 2) of the experimental data array. 
measuredSpecies = {{'protein deGFP*'}};
msIx = 2; % this is the index of the measured species in the data array 
% from data_dsg2014. There are two species: 1: mRNA and 2: GFP. 


%% setup the MCMC simulation parameters
stdev = 1; % i have no idea what a good value is
tightening = 1; % i have no idea what a good value is
nW = 400; % actual: 200 - 600 ish
stepsize = 1.3; % actual: 1.1 to 4 ish
niter = 10; % actual: 2 - 30 ish,
npoints = 4e4; % actual: 2e4 to 2e5 ish (or even 1e6 of the number of 
%                        params is small)
thinning = 4; % actual: 10 to 40 ish

%% pull all this together into an output struct. 
% the mcmc info struct now is an array struct, the way struct should be used!

runsim_info = struct('stdev', {stdev}, ...
    'tightening', {tightening}, ...
    'nW', {nW}, ...
    'stepSize', {stepsize}, ...
    'nIter', {niter}, ...
    'nPoints', {npoints}, ...
    'thinning', {thinning}, ...
    'parallel', true);

% for now we simply make the model_info have just one model (topology). 
% But the code will be written in a way such that multiple models can be used. 

model_info = struct(...
    'circuitInfo',{circuitInfo},...
    'modelObj', {modelObj},... % array of model objects (different topologies)
    'modelName', {modelObj.name},...; % model names. 
    'namesUnord', {activeNames}, ... % names of parameters per model, unordered. 
    'paramMaps', {paramMap}, ... % each paramMap is a matrix mapping models to the master vector. 
    'dosedNames', {dosedNames},... % cell arrays of species. cell array corresponds
     ...                               % to a model.
    'dosedVals', {dosedVals},...  % matrices of dose vals 
    'measuredSpecies', {measuredSpecies}, ... % cell array of cell arrays of 
                      ...                  % species names. the elements of the inner
                      ...                  % cell array get summed. 
    'measuredSpeciesIndex', {msIx},...                                    
    'dataToMapTo', dataIndices); % each dataToMapTo property within an element of the 
                            % model_info array is a vector of length # of geometries. 
                            
                            
semanticGroups = num2cell((1:10)'); %arrayfun(@num2str, 1:10, 'UniformOutput', false);


estParamsIx = setdiff((1:length(masterVector))', fixedParams);

%% master parameter vector, param ranges, 
master_info = struct(...
    'estNames', {estParams},...
    'masterVector', {masterVector},... % master vector is the set of all parameters that get distributed 
     ...                                   % across all models and topologies. not just the estimated params
    'paramRanges', {paramRanges},... % paramRanges are the ranges of parameter values for the values in the 
    ... % master vector that are estimated. ie, are specified by estParamIx as defined above. 
    'fixedParams', {fixedParams},...   % indexes of the parameters that are fixed (withing master vector)
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
                                        % then all parameters are assumed to be in their distinct
                                        % groups. 


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



mcmc_info = struct('runsim_info', runsim_info, ...
    'model_info', model_info,...
    'master_info', master_info);

end