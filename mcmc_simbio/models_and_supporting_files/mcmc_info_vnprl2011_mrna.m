function mcmc_info = mcmc_info_vnprl2011_mrna(modelObj)
    % version2 mcmc_info struct. Compatible with multimodel parameter inference. 
    % This mcmc info has two linked estimation problems: 
    % 1) transcription estimation 
    % 2) RNA degradation 
    % 
    % There is a second file in this series, mcmc_info_vnprl2011_protein,
    % that tries to fit the protein data, with the mRNA parameters fixed to
    % a few values found in this estimation. 
    % 
    % Finally there is a third file in this series that starts from all the
    % parameters estimated in both the first and second estimations, and
    % tries to fit all the parameters to all the data simultaneously within
    % a relatively narrow range around the pre-fit parameters. 
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
    circuitInfo1 = ...
    ['This is a simple constitutive gene expression model \n'...
    'built using the TXTL modeling toolbox. It models DNA binding \n'...
    'to RNAP and nucleotides, followed by transcription. The resulting\n'...
    'mRNA can degrade and participate in translation. The former is \n'...
    'modeled as a enzymatic reaction involving every complex containing \n'...
    'mRNA. The latter involves binding to Ribosomes, followed by amino acids \n'...
    'and ATP, and finally elongation and termination resulting in protein.'];

    circuitInfo2 = ...
    ['This is simple enzymatic rna degradation. Note that here mRNA can \n'...
    'to ribosomes and other species, but these offer no protection. '...
    'from rna degradation. \n '];

    

    %{
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
    %}
    
    % names of the parameters and species to set allow for setting in the
    % exported model 
    %%
    activeNames1 = {...
        'TX_elong_glob'                     1           [0.5 10]     
        'AGTPdeg_time'                      7200        [1800 18000]
        'AGTPdeg_rate'                      0.0002      [1e-5 1e-1]
        'TXTL_P70_RNAPbound_Kd'             12          [0.1 1000]
        'TXTL_P70_RNAPbound_F'              17          [0.1 100]
        'TXTL_RNAPBOUND_TERMINATION_RATE'   0.07        [1e-4 10]
        'TXTL_NTP_RNAP_1_Kd'                2           [0.1 1000]
        'TXTL_NTP_RNAP_1_F'                 10          [0.1 100]
        'TXTL_NTP_RNAP_2_Kd'                10          [0.1 1000]
        'TXTL_NTP_RNAP_2_F'                 1           [0.1 100]
        'TXTL_RNAdeg_Kd'                    2000        [1 5000]
        'TXTL_RNAdeg_F'                     1           [0.1 1000]
        'TXTL_RNAdeg_kc'                    0.001       [1e-4 1]
        'RNAP'                              30          [5 500]
        'RNase'                             200         [10 1000]}; 
    
    activeNames2 = {...
        'TXTL_RNAdeg_Kd'                    2000        [1 5000]
        'TXTL_RNAdeg_F'                     1           [0.1 1000]
        'TXTL_RNAdeg_kc'                    0.001       [1e-4 1]
        'RNase'                             200         [10 1000]}; 
    %%
    % Names of parameters and species to actually estimate. 
    estParams = activeNames1(:,1);   
    % fixedParams vector
    fixedParams = []; % none
    masterVector = zeros(length(activeNames1(:,1)), 1); % log transformed. 
    
    % paramMap is a matrix mapping the parameters in the master vector to the 
    % (unordered) list of parameters in the model. (obvioulsy within the code 
    % these parameters get ordered before they are used in the exported model)
    % More precisely, Let pp = paramMap(:, 1); then masterVector(pp) is the list 
    % of parameters for the first geometry within that topology. 
    % One such matrix exists for each topology. It has dimnesions 
    % length(model_info(i).namesUnord) x number of geometries associated with that topo. 
    paramMap1 = [1:length(activeNames1(:,1))]';
    paramMap2 = [11 12 13 15]';

    % parameter ranges (for the to-be-estimated parameters in the master
    % vector)
    paramRanges = log(cell2mat(activeNames1(:,3)));
    
    
%% next we define the dosing strategy. 
    dosedNames1 = {'DNA p70--utr1--deGFP'};
    dosedVals1 = [0.5 2 5 20];

    dosedNames2 = {'RNA utr1--deGFP'};
    dosedVals2 = [37.5 75 150 200 600 700 800 900 1000];
%% create the measured species cell array
measuredSpecies = {{'[RNA utr1--deGFP]',...
    '[Ribo:RNA utr1--deGFP]',...
    '[AA:2AGTP:Ribo:RNA utr1--deGFP]', ...
    '[term_Ribo:RNA utr1--deGFP]',...
    '[AA:Ribo:RNA utr1--deGFP]'...
    '[RNA utr1--deGFP:RNase]',...
    '[Ribo:RNA utr1--deGFP:RNase]',...
    '[AA:2AGTP:Ribo:RNA utr1--deGFP:RNase]', ...
    '[term_Ribo:RNA utr1--deGFP:RNase]',...
    '[AA:Ribo:RNA utr1--deGFP:RNase]'}};
msIx = 1; % this is the index of the measured species in the data array 
% from data_dsg2014. There are two species: 1: mRNA and 2: GFP. 


%% setup the MCMC simulation parameters
stdev = 1; % i have no idea what a good value is
tightening = 1; % i have no idea what a good value is
nW = 400; % actual: 200 - 600 ish
stepsize = 1.5; % actual: 1.1 to 4 ish
niter = 20; % actual: 2 - 30 ish,
npoints = 4e4; % actual: 2e4 to 2e5 ish (or even 1e6 of the number of 
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
    'parallel', true);

model_info = struct(...
    'circuitInfo',{circuitInfo1, circuitInfo2},...
    'modelObj', {modelObj,modelObj},... % array of model objects (different topologies)
    'modelName',  modelObj.name,...; % model names. 
    'namesUnord', {activeNames1(:,1),activeNames2(:,1)}, ... % names of parameters per model, unordered. 
    'paramMaps', {paramMap1, paramMap2}, ... % paramMap is a matrix mapping models to master vector. 
    'dosedNames', {dosedNames1, dosedNames2},... % cell arrays of species. cell array corresponds
     ...                               % to a model.
    'dosedVals', {dosedVals1, dosedVals2},...  % matrices of dose vals 
    'measuredSpecies', {measuredSpecies, measuredSpecies}, ... % cell array of cell arrays of 
                      ...                  % species names. the elements of the inner
                      ...                  % cell array get summed. 
    'measuredSpeciesIndex', msIx,...  % maps measuredSpecies to the species in data array
    'dataToMapTo', {1, 3}); % each dataToMapTo property within an element of the 
                            % model_info array is a vector of length # of geometries. 
    % data indices tell us which data set to use for each topology (model) - geometry pair
    % from the data_info struct array. 
    
                            
semanticGroups = num2cell((1:length(estParams))'); 
%arrayfun(@num2str, 1:10, 'UniformOutput', false);
estParamsIx = setdiff((1:length(masterVector))', fixedParams);

%% master parameter vector, param ranges, 
master_info = struct(...
    'estNames', {estParams},...         
    'masterVector', {masterVector},...  
    'paramRanges', {paramRanges},...    % 
    'fixedParams', {fixedParams},...    % indexes of the fixed params (withing master vector)
    'semanticGroups', {semanticGroups});% EITHER EMPTY OR
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

end