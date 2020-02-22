function di = data_ZSIFFL
% IFFL data collected by ZSun and SGuo
% 
% This data info is of length 5: i.e., there are 5 datasets in it. 
%     pTet constitutive data
%     TetR repression data
%     aTc induction data
%     pLac constitutive data
%     3oc12 AHL induction data
% 
% Copyright (c) 2018, Vipul Singhal, Caltech
% Permission is hereby granted, free of charge, to any person obtaining a copy
% of this software and associated documentation files (the "Software"), to deal
% in the Software without restriction, including without limitation the rights
% to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
% copies of the Software, and to permit persons to whom the Software is
% furnished to do so, subject to the following conditions:
%
% The above copyright notice and this permission notice shall be included in all
% copies or substantial portions of the Software.
%
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
% IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
% FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
% AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
% LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
% OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
% SOFTWARE.

dataPath = [pwd '/mcmc_simbio/exp_data/Zach_IFFL_raw'];
%% Step 1: Build time vectors and data matrices

tv =  60*(0:8:(8*60))';
dmexp = zeros(length(tv), 3*35);% data matrix of experimental data in triplicate

% load S1 data
load([dataPath '/S1/S1data.mat'], 'Set1Data');
for i = 1:7
    % trajectories of original data in triplicate
    [T, it_exp, ~] = unique(Set1Data(:,1,i));
    dmexp(:,3*(i-1)+(1:3))= interp1(T, Set1Data(it_exp,2:4,i),tv);
end

% load S2 data 
load([dataPath '/S2/S2data.mat'], 'Set2Data');
for i = 1:7
    offset = 7;
    [T, it_exp, ~] = unique(Set2Data(:,1,i));
    dmexp(:,3*(offset+(i-1))+(1:3))= interp1(T, Set2Data(it_exp,2:4,i),tv);
end

% load S3 and S4 data 
load([dataPath '/S3/S3S4data.mat'], 'Set3Data', 'transformedSet4data');

for i = 1:7
    offset = 14;
    [T, it_exp, ~] = unique(Set3Data(:,1,i));
    dmexp(:,3*(offset+(i-1))+(1:3))= interp1(T, Set3Data(it_exp,2:4,i),tv);
    
end

for i = 1:7
    offset = 21;
    [T, it_exp, ~] = unique(transformedSet4data(:,1,i));
    dmexp(:,3*(offset+(i-1))+(1:3))= interp1(T, transformedSet4data(it_exp,2:4,i),tv);
end

% load S5 data 
%load('/Users/vipulsinghal/Dropbox/110114/Vipul_computational/Computational Tools/S5/S5data.mat');
load([dataPath '/S5/S5data.mat'], 'Set5Data');
for i = 1:7
    offset = 28;
    % trajectories of original data in triplicate
    [T, it_exp, ~] = unique(Set5Data(:,1,i));
    dmexp(:,3*(offset+(i-1))+(1:3))= interp1(T, Set5Data(it_exp,2:4,i),tv);
end
%% Step 2: build index arrays

% experimental data
ixexp = {[1 2 3; 4 5 6; 7 8 9; 10 11 12; 13 14 15; 16 17 18; 19 20 21], []; ... % 7 conditions for the ptet constitutive. 3 repeats each. 
     7*3+[1 2 3; 4 5 6; 7 8 9; 10 11 12; 13 14 15; 16 17 18; 19 20 21], []; 
    14*3+[1 2 3; 4 5 6; 7 8 9; 10 11 12; 13 14 15; 16 17 18; 19 20 21], []; ... % 7 conditions too... and so on. 
    21*3+[1 2 3; 4 5 6; 7 8 9; 10 11 12; 13 14 15; 16 17 18; 19 20 21], [];
    28*3+[1 2 3; 4 5 6; 7 8 9; 10 11 12; 13 14 15; 16 17 18; 19 20 21], []};

%% Step 3: build title arrays
datadescription = {'Zach IFFL training data: pTet constitutive'; 
    'Zach IFFL training data: pLac constitutive';
    'Zach IFFL training data: TetR repression'; 
    'Zach IFFL training data: aTc induction';
    'Zach IFFL training data: 3OC12 AHL induction'};
% final order is time x species x replicates x doses
% da = permute(da, [1, 3, 4, 2]); 
dimlabels = {'time points', 'measured species', 'replicates', 'doses'};

%%
dv1 = [4        2       1       .5      .25     .125    .0625];
dv2 = [2        1       .5      .25     .125    .0625   .03125];
dv3 = [2        .2      .02     .002    .0002   .00002  0];
dv4 = [10000    1000    100     10      1       .1      0];
dv5 = [10000    1000    100     10      1       .1      0];

mn = {'deGFP'};
dosedNames1 = {'pTet-deGFP DNA'};
dosedNames2 = {'pLac-deGFP DNA'};
dosedNames3 = {'pLac-TetR DNA'};
dosedNames4 = {'aTc'};
dosedNames5 = {'3OC12HSL'};
% 
da = cell(5, 1);
for i = 1:5
    
    da{i} = zeros(length(tv), 1, 3, 7);
    for rep = 1:3
        for doseid = 1:7
            da{i}(:, 1, rep, doseid) = 1000*dmexp(:, ixexp{i,1}(doseid, rep));
        end
    end
    
    
end

       

        di1 = struct('dataInfo', datadescription(1), ...
			'timeVector', {tv}, ...
			'timeUnits', {'seconds'},...
			'dataArray', da(1),...
			'measuredNames', {mn},...
			'dataUnits', {{'nM'}},...
			'dimensionLabels', {dimlabels}, ...
			'dosedNames', {dosedNames1},...
			'dosedVals', {dv1}, ...
			'doseUnits', 'nM');
        
        di2 = struct('dataInfo', datadescription(2), ...
			'timeVector', {tv}, ...
			'timeUnits', {'seconds'},...
			'dataArray', da(2),...
			'measuredNames', {mn},...
			'dataUnits', {{'nM'}},...
			'dimensionLabels', {dimlabels}, ...
			'dosedNames', {dosedNames2},...
			'dosedVals', {dv2}, ...
			'doseUnits', 'nM');
        di3 = struct('dataInfo', datadescription(3), ...
			'timeVector', {tv}, ...
			'timeUnits', {'seconds'},...
			'dataArray', da(3),...
			'measuredNames', {mn},...
			'dataUnits', {{'nM'}},...
			'dimensionLabels', {dimlabels}, ...
			'dosedNames', {dosedNames3},...
			'dosedVals', {dv3}, ...
			'doseUnits', 'nM');
        di4 = struct('dataInfo', datadescription(4), ...
			'timeVector', {tv}, ...
			'timeUnits', {'seconds'},...
			'dataArray', da(4),...
			'measuredNames', {mn},...
			'dataUnits', {{'nM'}},...
			'dimensionLabels', {dimlabels}, ...
			'dosedNames', {dosedNames4},...
			'dosedVals', {dv4}, ...
			'doseUnits', 'nM');
        di5 = struct('dataInfo', datadescription(5), ...
			'timeVector', {tv}, ...
			'timeUnits', {'seconds'},...
			'dataArray', da(5),...
			'measuredNames', {mn},...
			'dataUnits', {{'nM'}},...
			'dimensionLabels', {dimlabels}, ...
			'dosedNames', {dosedNames5},...
			'dosedVals', {dv5}, ...
			'doseUnits', 'nM');

  
 di = [di1, di2, di3, di4, di5];
 
end





