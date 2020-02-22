function mcmc_info = mcmc_info_test015_corr1_Ffix(modelObj, rkcp, cpol)
% Single topology: tetR repression circuit
% 2 geometries (extracts)
% kf(s) fixed. 
% Separate estimation of ESPs and CSPs. 
% 
% ~~~ MODEL ~~~
% D_T + P <-> D_T:P -> D_T + P + T
% D_G + P <-> D_G:P -> D_G + P + G
% 2 T <-> T2
% D_G + T2 <-> D_G:T2
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
    [  'D_T + P <-> D_T:P -> D_T + P + T (kfdT, krdT, kcpT) \n'...
       'D_G + P <-> D_G:P -> D_G + P + G (kfdG, krdG, kcpG) \n'...
       '2 T <-> T2 (kfdim, krdim)\n'...
       'D_G + T2 <-> D_G:T2 (kfseq, krseq) \n'...
    'single topology, single geometry.'];

% good starting values? seems like it!
%    SimBiology Parameter Array
% 
%    Index:    Name:       Value:    ValueUnits:
%    1         kfdT        1         
%    2         krdT        6         
%    3         kcp         0.012     
%    4         kfdG        1         
%    5         krdG        6         
%    6         kfdimTet    2         
%    7         krdimTet    4         
%    8         kfseqTet    2         
%    9         krseqTet    4         
% 
% mobj.Parameters.Value
% Undefined function or variable 'Value'.
%  
% mobj.Parameters(1).Value
% 
% ans =
% 
%      1
% 
% log([1 6 .012 1 6 2 4 2 4 100])
% 
% ans =
% 
%          0    1.7918   -4.4228         0    1.7918    0.6931    1.3863    0.6931    1.3863    4.6052
         
%cpol = exp(1.0212); % nM
rkfdG = 1; % nM-1s-1
rkrdG = 6; % s-1
rkfdT = 1; % 
rkrdT = 6;
%rkcp = exp( -7.6000); %s-1
rkfdimTet = 2; % nM-1s-1 
rkrdimTet = 4; % s-1
rkfseqTet = 2; % nM-1s-1 
rkrseqTet = 4; % s-1

activeNames = ...
    {'kfdG'
    'krdG'
    'kfdT'
    'krdT'
    'kfdimTet'
    'krdimTet'
    'kfseqTet'
    'krseqTet'
    'kcp'
    'pol'};

masterVector = log([...
rkfdG 
rkrdG
rkfdT
rkrdT
rkfdimTet
rkrdimTet
rkfseqTet
rkrseqTet
rkcp
cpol]);
% fixedParams vector
fixedParams = [5 9 10];

estParamsIx = setdiff((1:length(masterVector))', fixedParams);

estParams = activeNames(estParamsIx); %{'kfdG'
%     'krdG'
%     'krdT'
%     'krdimTet'
%     'krseqTet'};

paramMap = [1:length(masterVector)]';
paramRanges =  [masterVector(estParamsIx)-10 masterVector(estParamsIx)+5];


dataIndices = [6];

%% next we define the dosing strategy.
dosedNames = {'dT'; 'dG'};
dosedVals = flipud([   10.0000   10.0000   10.0000   10.0000   10.0000   10.0000   10.0000   10.0000
         0    0.2500    0.5000    0.7500    1.0000    2.0000    5.0000   10.0000]);

         
measuredSpecies = {{'pG'}};
msIx = 1; %

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