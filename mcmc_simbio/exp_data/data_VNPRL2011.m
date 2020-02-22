function di = data_VNPRL2011
% This data is inspired from VNPRL 2011: 
% Karzbrun, Eyal, Jonghyeon Shin, Roy H. Bar-Ziv, and Vincent Noireaux. 
% ?Coarse-Grained Dynamics of Protein Synthesis in a Cell-Free System.? 
% Physical Review Letters 106, no. 4 (January 24, 2011): 48104. 
% https://doi.org/10.1103/PhysRevLett.106.048104.
% 
% It is not exactly the data from this paper, instead is created
% artificially to reflect the conclusions of that work. 
% Some of the main conclusions of that paper were: 
% - Transcriptional Elongation rate: 1 ntp/s
% - Translational Elongation rate: >4 aa/s
% - mRNA exponential decay, even when purified RNA is in excess of 200nM
% - mRNA degradation half life: 10 - 14 min
% - 30nM RNAP conc
% - 1.5nM RNAP - promoter Kd
% - Protein production linear in mRNA (TL machinery not saturated)
% - 1uM protein in 1h, And anywhere between 3 to 10 uM by the end (5 ish
% hours)
% - dp/dt|max is about 30 to 40 nM / min for proteins
% - For 30 nM of DNA, mRNA steady state is 20 - 30 nM. 
% 
% The estimation uses the following fits: 
% - mrna data from ACS DSG with the conc values all divided by 10 (the
% original values in that paper seem unlikey, since they would use too many
% nucleotides. 
% - protein max values come down from 10uM to 10uM
% - the rna degradation values just basically are exactly the same as in 
% the ACS paper 

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



%
% define the doses: DNA at 0.5, 2, 5 and 20nM
dosedNames = {'GFP DNA'};
dv = [0.5, 2, 5, 20]; % dose vals of the dna in nM

% get the data and the time vector (in seconds)
[tmg, ymg, metamg] = load_ACSDSG2014('MGapt'); % y is ntimepoints x ndoses
[tgfp, ygfp, metagfp] = load_ACSDSG2014('deGFP');
[trnadeg, yrnadeg, metarnadeg] = load_ACSDSG2014('RNAdeg');
[tgfp_deriv, ygfp_deriv, metagfp_deriv] = load_ACSDSG2014('deGFP_deriv');
ygfp = ygfp(:,[2 4 5 6]);

tv = tmg;
da = zeros(length(tv), 4, 2, 1);
da(:, :, 1, 1) = ymg/10; % the mg aptameter data
da(:, :, 2, 1) = ygfp/1.8; % the gfp data

mn = {'MG aptamer', 'deGFP'};
% final order is time x species x replicates x doses
da = permute(da, [1, 3, 4, 2]); 

dimlabels = {'time points', 'measured species', 'replicates', 'doses'};
datadescription = ...
['Data from ACSDSG 2014, modified according to VNPRL 2011, \n '...
'MG aptamer as the first measured species, \n' ...
'GFP as the second measured species. Dosing at 0.5, 2, 5, 20nM'...
];

di1 = struct('dataInfo', {datadescription}, ...
			'timeVector', {tv}, ...
			'timeUnits', {'seconds'},...
			'dataArray', {da},...
			'measuredNames', {mn},...
			'dataUnits', {{'nM', 'nM'}},...
			'dimensionLabels', {dimlabels}, ...
			'dosedNames', {dosedNames},...
			'dosedVals', {dv}, ...
			'doseUnits', 'nM');
        
ygfp_deriv = ygfp_deriv(:,[2 4 5 6]);

da = zeros(length(tgfp_deriv), 4, 1, 1);
da(:, :, 1, 1) = ygfp_deriv/1.8; % the gfp data
da = permute(da, [1, 3, 4, 2]);

datadescription = ...
['Data from ACSDSG 2014, modified according to VNPRL 2011, \n '...
'dGFP/dt as the measured species. Dosing at 0.5, 2, 5, 20nM'...
];

di2 = struct('dataInfo', {datadescription}, ...
			'timeVector', {tgfp_deriv}, ...
			'timeUnits', {'seconds'},...
			'dataArray', {da},...
			'measuredNames', {mn},...
			'dataUnits', {{'nM'}},...
			'dimensionLabels', {dimlabels}, ...
			'dosedNames', {dosedNames},...
			'dosedVals', {dv}, ...
			'doseUnits', 'nM');           
        
 mn = {'MG aptamer'};
dv = [37.5 75 150 200 600 700 800 900 1000];
dosedNames = {'Purified RNA'};
da = zeros(length(trnadeg), 9, 1, 1);
da = yrnadeg;
da = permute(da, [1, 3, 4, 2]); 
datadescription = ...
['Data from ACSDSG 2014, modified according to VNPRL 2011, \n '...
'MG aptamer as the measured species, \n' ...
'Purified mRNA added. Degradation of mRNA observed. '...
'Dosing: 37.5 75 150 200 600 700 800 900 1000 nM of mRNA'...
];

di3 = struct('dataInfo', {datadescription}, ...
			'timeVector', {trnadeg}, ...
			'timeUnits', {'seconds'},...
			'dataArray', {da},...
			'measuredNames', {mn},...
			'dataUnits', {{'nM'}},...
			'dimensionLabels', {dimlabels}, ...
			'dosedNames', {dosedNames},...
			'dosedVals', {dv}, ...
			'doseUnits', 'nM'); 
        
 di = [di1, di2, di3];
 
end





