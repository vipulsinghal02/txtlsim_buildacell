function mcmc_info = mcmc_info_constgfp5ii(modelObj)
% Single topology: constitutive gfp circuit
% 2 geometries (extracts)
% kf(s) fixed. 
% separate ESPs and joint CSPs. 
%
% MODEL: tetRrepression 2.
% dG + pol <-> dG_pol -> dG + pol + mG
% mG + ribo <-> mG_ribo -> mG + ribo + pG
% mG -> null
%

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
    ['dG + pol <-> dG_pol -> dG + pol + mG  \n'...
    'mG + ribo <-> mG_ribo -> mG + ribo + pG  \n'...
    'mG -> null  \n'];



rkfdG = 1; % nM-1s-1
rkrdG = 60; % s-1

rkfpG = 2; % nM-1s-1
rkrpG = 60; % s-1

cpol1 = 100; % nM
cribo1 = 50; %nM
rkcm1 = 0.001; %s-1
rkcp1 = 1/36;
rdel_m1 = log(2)/720; % 12 min half life of mrna

cpol2 = cpol1*2;
cribo2 = cribo1*2;
rkcm2 = rkcm1*2;
rkcp2 = rkcp1*2;
rdel_m2 = rdel_m1*2;


activeNames = ...
    {'kfdG'
    'krdG'
    'kfpG'
    'krpG'
    'kcm'
    'kcp'
    'del_m'
    'pol'
    'ribo'};

masterVector = log([...
    rkfdG
    rkrdG
    rkfpG
    rkrpG
    rkcm1
    rkcp1
    rdel_m1
    cpol1
    cribo1
    rkcm2
    rkcp2
    rdel_m2
    cpol2
    cribo2]);
% fixedParams vector
fixedParams = [1 3];

estParamsIx = setdiff((1:length(masterVector))', fixedParams);

estParams = {'krdG'
    'krpG'
    'kcm1'
    'kcp1'
    'del_m1'
    'pol1'
    'ribo1'
    'kcm2'
    'kcp2'
    'del_m2'
    'pol2'
    'ribo2'};

paramMap1 = [1:4 5:9]';
paramMap2 = [1:4 10:14]';
paramMap = [paramMap1, paramMap2];
paramRanges =  [masterVector(estParamsIx)-5 masterVector(estParamsIx)+3];


dataIndices = [1 2];

%% next we define the dosing strategy.
dosedNames = {'dG'};
dosedVals = [10 30 60];
measuredSpecies = {{'pG'}};
msIx = 1; %

%% setup the MCMC simulation parameters
stdev = 1; % i have no idea what a good value is
tightening = 1; % i have no idea what a good value is
nW = 400; % actual: 200 - 600 ish
stepsize = 1.1; % actual: 1.1 to 4 ish
niter = 40; % actual: 2 - 30 ish,
npoints = 8e4; % actual: 1e4 to 2e5 ish (or even 1e6 of the number of
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