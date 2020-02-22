function di = data_VS2017_deGFP_MGa
% This dataset was collected by Vipul Singhal on 07-03-2017 (MM-DD-YYYY)
% and is stored in 
% /Users/vipulsinghal/Dropbox/Documents/a_Research/Labwork
% /Lab_2017/2017_07_13_prom_saturation_test
% It shows deGFP-MGa curves in 4 different extracts. 
% This dataset normalizes the RNA levels using the ACS DGG paper as
% follows: Use the biotek calibrations to get the nM values for 
% the protein
% conc in each extract. Then use the ACS DSG data to set the MG apramer
% peak levels. for example, if we have ( 20uM, 400nM) from the ACS DSG
% data, and (10uM, 10k AFU) from the eVS data, then the eVS MGA number is
% 10/20 * 400 = 200nM. and all the mRNA alibrations are going to follow
% this: 200/10000 = 0.02nM / AFU mRNA calibration. Can check if this
% matches
% for all the concentrations and extracts. (or at least matches within an
% extract). This scheme is basically assuming that at 20nM of DNA, the
% extract used by DSG and by VS have the same mapping from mRNA to GFP
% (even if the total levels produced are different)
% 
% Copyright (c) 2018, Vipul Singhal, Caltech
% Permission is hereby granted, free of charge, to any person obtaining
% a copy
% of this software and associated documentation files (the "Software"),
% to deal
% in the Software without restriction, including without limitation 
% the rights
% to use, copy, modify, merge, publish, distribute, sublicense, 
%and/or sell
% copies of the Software, and to permit persons to whom the Software is
% furnished to do so, subject to the following conditions:
%
% The above copyright notice and this permission notice shall be 
%included
% in all
% copies or substantial portions of the Software.
%
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, 
%EXPRESS OR
% IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
%MERCHANTABILITY,
% FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. 
%IN NO EVENT SHALL THE
% AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, 
%DAMAGES OR OTHER
% LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, 
%ARISING FROM,
% OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER 
%DEALINGS IN THE
% SOFTWARE.


% start by gettin the raw data
addpath(['/Users/vipulsinghal/Dropbox/Documents'...
    '/a_Research/Labwork/Lab_2017/'...
    '2017_07_13_prom_saturation_test/'])
[dFull, t] = data170713_saturation;
dG = dFull{1}; % already in uM. 
dMGA=dFull{5};
% subtract negative controls
origIgfp = [19+[1:11 20:26] 19+38+[1:11 20:26]...
    19+2*38+[1:11 20:26] 19+3*38+[1:11 20:26]];
subIgfp = [ 19+19*ones(1,18) 19+38+19*ones(1,18)...
    19+2*38+19*ones(1,18) 19+3*38+19*ones(1,18)];

origImga = [19+[12:17] 19+38+[12:17] 19+2*38+[12:17] 19+2*38+[12:17]];
subImga = [[19+18*ones(1,6) 19+38+18*ones(1,6)...
    19+2*38+18*ones(1,6) 19+3*38+18*ones(1,6)]];

dG = subtractCols(dG, origIgfp, subIgfp, true);
dMGA = subtractCols(dMGA, origImga, subImga, true);

ixexp4 = {19+11+[1  ; 2  ; 3  ; 4 ; 5 ; 6 ],...
   19+11+37+1+    [1  ; 2  ; 3  ; 4 ; 5 ; 6 ], ...
   19+11+74+2+    [1  ; 2  ; 3  ; 4 ; 5 ; 6 ], ...
    19+11+111+3+  [1  ; 2  ; 3  ; 4 ; 5 ; 6 ]};

%%

% 4 different extracts
nExtracts = length(ixexp4);

% define the doses: DNA at 0.5, 2, 5 and 20nM
dosedNames = {'GFP MGa DNA'};
dv = [5, 10, 15, 20, 30, 40]; % dose vals of the dna in nM

% construct the data array objects for the four extracts. 
tv = t;
da = zeros(length(t), size(dv, 2), 2, 1); 
% to populate using the for loop, start with the order nTime x nDoses x nMS x nRep 
% THEN reorder using permute. 
dacell = cell(1, nExtracts);

% renormalization factors for MGa data. 
% The max expression of GFP at 20nM DNA in the 2014 ACS paper was 18079nM 
% and the corresponding max in the RNA levels was 388.2915 nM
% lets assume that eVS has the same characteristic at 20nM DNA. We will use
% this to get the platereader calibration for AFU -> nM. 

deGFP_eVS_20nM_max = max(dG(:,ixexp4{1}(4)));
deGFP_ACS_20nM_max = 18079;
MGa_ACS_20nM_max = 388.2915;
MGa_eVS_20nM_max = deGFP_eVS_20nM_max/deGFP_ACS_20nM_max*MGa_ACS_20nM_max;
MGa_eVS_20nM_AFU = max(dMGA(:,ixexp4{1}(4)));
MGA_nM_per_AFU = MGa_eVS_20nM_max/MGa_eVS_20nM_AFU;


for i = 1:nExtracts
    da(:,:,1,1) = dMGA(:,ixexp4{i})*MGA_nM_per_AFU;
    da(:,:,2,1) = dG(:,ixexp4{i});
    daord = permute(da, [1, 3, 4, 2]);
    dacell{i} = daord;
end

mn = {{'MG aptamer'}, {'deGFP'}};

dimlabels = {'time points', 'measured species', 'replicates', 'doses'};
datadescription1 = ...
['Data from experiments done by VS, July, 2017 \n '...
'EXTRACT: eVS \n '...
'MG aptamer as the first measured species, \n' ...
'GFP as the second measured species. Dosing at 5, 10, 15, 20, 30, 40nM '...
];
datadescription2 = ...
['Data from experiments done by VS, July, 2017 \n '...
'EXTRACT: eSG \n '...
'MG aptamer as the first measured species, \n' ...
'GFP as the second measured species. Dosing at 5, 10, 15, 20, 30, 40nM '...
];
datadescription3 = ...
['Data from experiments done by VS, July, 2017 \n '...
'EXTRACT: eJP \n '...
'MG aptamer as the first measured species, \n' ...
'GFP as the second measured species. Dosing at 5, 10, 15, 20, 30, 40nM '...
];
datadescription4 = ...
['Data from experiments done by VS, July, 2017 \n '...
'EXTRACT: eGA \n '...
'MG aptamer as the first measured species, \n' ...
'GFP as the second measured species. Dosing at 5, 10, 15, 20, 30, 40nM '...
];

di = ...
    struct(...
'dataInfo', {datadescription1, datadescription2,datadescription3,datadescription4}, ...
'timeVector', {tv}, ...
'timeUnits', {'seconds'},...
'dataArray', dacell,...
'measuredNames', {mn},...
'dataUnits', {{'nM', 'nM'}},...
'dimensionLabels', {dimlabels}, ...
'dosedNames', {dosedNames},...
'dosedVals', {dv}, ...
'doseUnits', {'nM'});
       
%   mcmc_trajectories([], di, [], [], [], [], 'just_data_info', true);
  
  
end





