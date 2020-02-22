function di = data_dsg2014
% Generate the data_info struct containing the data from Figure 1 of 
% the 2014 ACS Synthetic Biology paper titled: 
% Gene Circuit Performance Characterization and Resource Usage in a 
% Cell-Free ‚ÄúBreadboard‚Äù by Siegal-Gaskins et al. This data involves a
% measurement of malachite green aptamer and green fluorescent protein 
% over a period of 800 minutes in the TX-TL cell free expression system. 
% The DNA initial conditions used are 0.5nM, 2nM, 5nM, and 20nM. 
% 
% ------------------------------------------
% The data_info struct has the fields: 
%
% 'dataInfo': A human readable text description of the data. 
% 
% 'timeVector': vector of timepoints
% 
% 'timeUnits': units of the time Vector
% 
% 'dataArray': An array contianing the raw data
% 
% 'measuredNames': A 1 x number of measured species cell array of the strings specifying 
% which species are dosed. These are not strings corresponding to the species 
% in the model. See mcmc_info constructor functions for that. 
% 
% 'dataUnits': A 1 x number measured species cell array of units corresponding to 
% the raw data in the dataArray
%
% 'dimensionLabels': a 1 by length(size(data_info.dataArray)) cell array of 
% labels for the dimensions of the dataArray. 
%
% 'dosedNames': A 1 x number of dosed species cell array of the strings specifying 
% which species are dosed. These are not strings corresponding to the species 
% in the model. See mcmc_info constructor functions for that. 
% 
% 'dosedVals': A matrix of dose values of size
% # of dosed species by # of dose combinations
% 
% 'doseUnits': A 1 x number of dosed species cell array of strings specifying the 
% units of the dosed species. 
%
% Note that data_info can also be an array of these inputs. This is useful in the 
% multi-modal mode (version 2 of the code).
% 
% 
% 
% ------------------------------------------

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

datadescription = ...
['Data from Figure 1 of the 2014 ACS Synthetic Biology paper titled:  \n '...
'Gene Circuit Performance Characterization and Resource Usage in a  \n '...
'Cell-FreeÄúBreadboardÄù by Siegal-Gaskins et al. This data involves a \n '...
'measurement of malachite green aptamer and green fluorescent protein  \n '...
'over a period of 800 minutes in the TX-TL cell free expression system.  \n '...
'The DNA initial conditions used are 0.5nM, 2nM, 5nM, and 20nM.' ];

% define the doses: DNA at 0.5, 2, 5 and 20nM
dosedNames = {'GFP DNA'};
dv = [0.5, 2, 5, 20]; % dose vals of the dna in nM

% get the data and the time vector (in seconds)
[tmg, ymg, ~] = load_ACSDSG2014('MGapt'); % y is ntimepoints x ndoses
[~, ygfp, ~] = load_ACSDSG2014('deGFP');
ygfp = ygfp(:,[2 4 5 6]);
tv = tmg*60;
da = zeros(length(tv), 4, 2, 1);
da(:, :, 1, 1) = ymg; % the mg aptameter data
da(:, :, 2, 1) = ygfp; % the gfp data

mn = {'MG aptamer', 'deGFP'};
% final order is time x species x replicates x doses
da = permute(da, [1, 3, 4, 2]); 

dimlabels = {'time points', 'measured species', 'replicates', 'doses'};

di = struct('dataInfo', {datadescription}, ...
			'timeVector', {tv}, ...
			'timeUnits', {'seconds'},...
			'dataArray', {da},...
			'measuredNames', {mn},...
			'dataUnits', {{'nM', 'nM'}},...
			'dimensionLabels', {dimlabels}, ...
			'dosedNames', {dosedNames},...
			'dosedVals', {dv}, ...
			'doseUnits', 'nM');

end





