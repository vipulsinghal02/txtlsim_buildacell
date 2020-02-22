function mcmc_info = mcmc_info_tierra2018_corr(modelObj, esps, csps)
% mcmc_info_tierra2018_corr.m 
% Correction: This is just an aTc induction model whose ESPs have been fixed by
% the calibration step and only the CSPs are estimated. 
% The model used is found in the .m file: model_aTc_induc1.m 
% There is no parameter sharing of any kind. 
% The data corresponds to the extract eEC5, and the experimental conditions 
% are 
% 4 nM ptetO1-deGFP       
% 0.4 nM placO1-tetR      
% 10, 1, 0.1, 0.01, 0.001 uM aTc       
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
    ['D_T + P <-> D_T:P -> D_T + P + T\n'...
     'D_G + P <-> D_G:P -> D_G + P + G\n'...
     '2 T <-> T2\n'...
     'D_G + T2 <-> D_G:T2\n'...
     '2 aTc <-> aTc2\n'...
     'aTc2 + T2 <-> aTc2:T2\n'...
     'single topology, single geometry.'];

cpol = esps(2); % nM
rkfdG = 0.5; % nM-1s-1
rkrdG = csps(1); % s-1
rkfdT = .2;
rkrdT = 30;
rkcp = esps(1); %s-1
rkfdimTet = .5; % nM-1s-1
rkrdimTet = 20; % s-1
rkfseqTet = .5; % nM-1s-1
rkrseqTet = 20; % s-1
% rkfdimaTc = .5;
% rkrdimaTc = 2;
rkfseqaTc = .5;
rkrseqaTc = 20;
activeNames = ...
    {'kfdG'
    'krdG'
    'kfdT'
    'krdT'
    'kfdimTet'
    'krdimTet'
    'kfseqTet'
    'krseqTet'
    'kfseqaTc'
    'krseqaTc'
    'kcp'
    'pol'};

%     'kfdimaTc'
%     'krdimaTc'
masterVector = log([...
rkfdG 
rkrdG
rkfdT
rkrdT
rkfdimTet
rkrdimTet
rkfseqTet
rkrseqTet
rkfseqaTc
rkrseqaTc
rkcp
cpol]);
% rkfdimaTc
% rkrdimaTc

% fixedParams vector
fixedParams = [1 3 5 7 9 11 12];

estParamsIx = setdiff((1:length(masterVector))', fixedParams);

estParams = {'krdG'
    'krdT'
    'krdimTet'
    'krseqTet'
    'rkrseqaTc'}; 
    % an 8 dimensional space gets searched. Initialize carefully (even 8 *can* be bad)
% 'rkrdimaTc'


paramMap = [1:length(masterVector)]';
paramRanges =  [masterVector(estParamsIx)-5 masterVector(estParamsIx)+5];


dataIndices = [1];

%% next we define the dosing strategy.
dosedNames = {'dT'; 'dG'; 'aTc2'};
dosedVals = [0.4 0.4 0.4 0.4 0.4; 
             4 4 4 4 4;
             10000 1000 100 10 1];

measuredSpecies = {{'pG'}};
msIx = 1; %

%% setup the MCMC simulation parameters
stdev = 1; % i have no idea what a good value is
tightening = .1; % i have no idea what a good value is
nW = 200; % actual: 200 - 600 ish
stepsize = 1.4; % actual: 1.1 to 4 ish
niter = 10; % actual: 2 - 30 ish,
npoints = 4e4; % actual: 2e4 to 2e5 ish (or even 1e6 of the number of
%                        params is small)
thinning = 20; % actual: 10 to 40 ish

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
    'paramMaps', {paramMap}, ... %paramMap: matrix mapping models to master vector.
    'dosedNames', {dosedNames},... % cell arrays of species. cell array corresponds
    ...                               % to a model.
    'dosedVals', {dosedVals},...  % matrices of dose vals
    'measuredSpecies', {measuredSpecies}, ... % cell array of cell arrays of
    ...                  % species names. the elements of the inner
    ...                  % cell array get summed.
    'measuredSpeciesIndex', {msIx},...
    'dataToMapTo', dataIndices); % each dataToMapTo property within an 
% element of the model_info array is a vector of length # of geometries.


semanticGroups = num2cell((1:length(estParams))'); 
%arrayfun(@num2str, 1:10, 'UniformOutput', false);


estParamsIx = setdiff((1:length(masterVector))', fixedParams);

%% master parameter vector, param ranges,
master_info = struct(...
    'estNames', {estParams},...
    'masterVector', {masterVector},...
    'paramRanges', {paramRanges},...
    'fixedParams', {fixedParams},...
    'semanticGroups', {semanticGroups});


mcmc_info = struct('runsim_info', runsim_info, ...
    'model_info', model_info,...
    'master_info', master_info);

end

