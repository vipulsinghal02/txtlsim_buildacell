function mcmc_info = mcmc_info_dsg2014_protein(modelObj)

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
'mRNA. The latter involves binging to Ribosomes, followed by amino acids \n'...
'and ATP, and finally elongation and termination resulting in protein.']

pcells = [{'TL_elong_glob'           4 
'TXTL_PROT_deGFP_MATURATION'         0.00231049 
'TXTL_UTR_UTR1_Kd'                   20     
'TXTL_UTR_UTR1_F'                    0.2 
'TL_AA_Kd'                           0.1           
'TL_AA_F'                            1             
'TL_AGTP_Kd'                         10            
'TL_AGTP_F'                          1             
'TXTL_RIBOBOUND_TERMINATION_RATE'    0.1}]; 
estNames = [pcells(:,1)
    {'Ribo'}]; % 3 8
estVals = [cell2mat(pcells(:,2))
    100];


% b = [estVals/10 estVals*10]  ;
% paramranges = log(b);

paramranges = [0 6
    -7 -3
    -1 7
    -3 3
    -1 6
    -1 6
    -1 6
    -1 6
    -6 1
    2 7];

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

%% setup the MCMC simulation parameters
stdev = 1; % i have no idea what a good value is
tightening = 1; % i have no idea what a good value is
nW = 100; % actual: 200 - 600 ish
stepsize = 3; % actual: 1.1 to 4 ish
niter = 2; % actual: 2 - 30 ish,
npoints = 3e3; % actual: 2e4 to 2e5 ish (or even 1e6 of the number of params is small)
thinning = 2; % actual: 10 to 40 ish

%% pull all this together into an output struct. 
mcmc_info = struct(...
    'circuitInfo',{circuitInfo},...
    'modelObj', {modelObj},...
    'modelname', {m.name},...;
    'namesUnord', {estNames}, ...
    'paramRanges', {paramranges},...
    'dosedNames', {dosedNames},...
    'dosedVals', {dosedVals},...
    'measuredSpecies', {measuredSpecies}, ...
    'stdev', {stdev}, ...
    'tightening', {tightening}, ...
    'nW', {nW}, ...
    'stepSize', {stepsize}, ...
    'nIter', {niter}, ...
    'nPoints', {npoints}, ...
    'thinning', {thinning}, ...
    'parallel', false);




end