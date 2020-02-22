function mcmc_info = mcmc_info_dsg2014_mrna_v2(modelObj)
    % Define the mcmc_info struct for the dsg2014 dataset, for the estimation 
    % of the mrna parameters, using only measurements of the mrna for 
    % the estimation. The data is from figure 1 of the paper:
    % Gene Circuit Performance Characterization and Resource Usage in a 
    % Cell-Free “Breadboard” by Siegal-Gaskins et al. 
    %
    % INPUTS: A Simbiology model object. 
    % 
    % OUTPUTS: You should set up an mcmc_info struct with the fields:
    %
    % The struct has fields: 
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

%% make a list of parameters and species to estimate 

% parameters and typical values. 
pcells = [{'TX_elong_glob',             5.1288      }   %   1 3
    {'AGTPdeg_time'                     7200        }   %   4 10
    {'AGTPdeg_rate'                     0.00022425  }   %   -10 -7
    {'TXTL_P70_RNAPbound_Kd'            12.39       }   %   -1 7
    {'TXTL_P70_RNAPbound_F'             17.001      }   %   1 7
    {'TXTL_RNAPBOUND_TERMINATION_RATE'  0.076026    }   %   -5 -1
    {'TXTL_NTP_RNAP_1_Kd'               2.142       }   %   -2 2
    {'TXTL_NTP_RNAP_1_F'                9.8854      }   %   1 4
    {'TXTL_NTP_RNAP_2_Kd'               13.33       }   %   1 4
    {'TXTL_NTP_RNAP_2_F'                0.75935     }   %   -2 2
    {'TXTL_RNAdeg_Kd'                   2185.8      }   %   5 10
    {'TXTL_RNAdeg_F'                      1         }   %   -5 2
    {'TXTL_RNAdeg_kc'                   0.0022382   }]; %   -7 -3

% species and typical values. 
estNames = [pcells(:,1)
    {'RNAP';... 
    'RNase'}]; 
estVals = [cell2mat(pcells(:,2))
    42.247
    81.795];

% parameter and species ranges to search, in log base e space. 
paramranges = [1 3
    4 10
    -10 -7
    -1 7
    1 7
    -5 -1
    -2 2
    1 4
    1 4
    -2 2
    5 10
    -5 2
    -7 -3
    1.5 6
    2 7];


%% Define the species to set the initial concentration of and what values to 
% vary it over. 
dosedNames = {'DNA p70--utr1--deGFP'};
dosedVals = [0.5 2 5 20];


%% create the measured species cell array
% this is a 1x2 cell array. each element of this cell array contains
% further cell arrays. The first such cell array is a list of all the bound
% and free versions of the RNA. The second cell array contains a single
% cell, which contains the GFP string. 
% all the species in the inner cell arrays get summed, and compared to the
% corresponding column (dimension 2) of the experimental data array. 
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

%% setup the MCMC simulation parameters
stdev = 1; % i have no idea what a good value is
tightening = 1; % i have no idea what a good value is
nW = 600; % actual: 200 - 600 ish
stepsize = 1.2; % actual: 2 to 4 ish
niter = 30; % actual: 2 - 20 ish,
npoints = 4e4; % actual: 1e5 ish
thinning = 10; % actual: 10 to 40 ish


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
    'parallel', true);






end