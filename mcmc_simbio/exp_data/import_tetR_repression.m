function [t , data_array, IC] = import_tetR_repression
% Function returns a data set as a 4d matlab array. 
% This data was presented at the 1st European Congress on Cell-Free 
% Synthetic Biology (ECCSB) held in Ascona, Switzerland in March 2017.
% It is stored in the folder  
% /Users/vipulsinghal/Dropbox/Documents/a_Research/Labwork/Lab_2017/2017_03_20
% The data describes an experiment where we vary the concentration of the
% plac-UTR1-tetR linear DNA at concentration values of, 0nM (unrepressed), 0.25, 0.5, 0.75,
% 1nM, 2nM, 5nM and 10nM, and measure the deGFP expression from a 5nM
% ptet-UTR1-deGFP linear DNA construct. 
% The deGFP fluorescence signal was measured on the Biotek 2 plate
% reader in the Murray lab at Caltech. 
%
% IC is a nDosedSpecies x nICs array
%
% the data is arranged as a matrix array as follows: 
% dim 1: time
% dim 2: dosing / ICs
% dim 3:  extract
% dim 4: species 

IC = [10 10 10 10 10 10 10 10; 
    0 0.25 0.5 0.75 1 2 5 10]; % first row is ptetGFP and second row is plac-tetR

% data_init ; % adds the directory with the function data_2017_03_20 in it. 

[dFull, t] = data_2017_03_20;
dG = dFull{1};
ntime = length(t); 
ndoses = 8;
nextracts = 3;
nspecies = 1;
% ptet GFP from 1 2 5 10 20 nM from G I K, 12 to 16, 
% 
% dvs2 ptetGFP 10nM, dvs3 tetR [.25 .5 .75 1 2 5 10]nM, IPTG 10uM, in G I K
% 17 to 20 and HJL 1 to 3. SKIPPED G21, so make sure that that column in
% removed from the data matrix. We also use the unrepressed index. 
sIx = alphanum2num('G1', 'L21', {'G12', 'I12', 'K12'; 'G17', 'I17', 'K17'; 'H1', 'J1', 'L1'; 'G14', 'I14', 'K14'});
eIx = alphanum2num('G1', 'L21', {'G16', 'I16', 'K16'; 'G20', 'I20', 'K20'; 'H3', 'J3', 'L3'});
d = dG;
lablFnt = 14;
titFnt = 16;
legFnt = 14;

% construct CUSTOMCOLIDX
CUSTOMCOLIDX = cell(2,3);
CUSTOMCOLIDX{1,1} = [sIx(1,1):eIx(1,1)];
CUSTOMCOLIDX{1,2} = [sIx(1,2):eIx(1,2)];
CUSTOMCOLIDX{1,3} = [sIx(1,3):eIx(1,3)];

CUSTOMCOLIDX{2,1} = [sIx(4,1) sIx(2,1):eIx(2,1) sIx(3,1):eIx(3,1)];
CUSTOMCOLIDX{2,2} = [sIx(4, 2) sIx(2,2):eIx(2,2) sIx(3,2):eIx(3,2)];
CUSTOMCOLIDX{2,3} = [sIx(4,3) sIx(2,3):eIx(2,3) sIx(3,3):eIx(3,3)];

% construct labels etc. 
titls = {'ptet-GFP varying, eVS', 'ptet-GFP varying, eMP', 'ptet-GFP varying, eSG';
   'ptet-GFP 5nM, pConst-tetR varying. eVS', 'ptet-GFP 5nM, pConst-tetR varying. eMP', 'ptet-GFP 5nM, pConst-tetR varying. eSG' };

legs = {{'ptetGFP 1nM', 'ptetGFP 2nM', 'ptetGFP 5nM', 'ptetGFP 10nM', 'ptetGFP 20nM'},...
    {'Un-repressed', 'pConst-tetR 0.25nM', 'pConst-tetR 0.5nM', 'pConst-tetR 0.75nM', 'pConst-tetR 1nM',...
    'pConst-tetR 2nM', 'pConst-tetR 5nM', 'pConst-tetR 10nM'}};
ylab = {'GFP, uM', 'GFP, uM'};
xlab = 'time, hours';
scrsz = get(groot,'ScreenSize');
figsz = [1 scrsz(4) scrsz(3) scrsz(4)]; %[left bottom width height].

% use the explicit version of the platereader plotting function. 
% 
% plot_platereader_subplots_explicit(dG(1:end-55,:),t(1:end-55,:)/3600, ...
%  CUSTOMCOLIDX, titls, legs, ylab,xlab, figsz, 14, 14, 16);
% 

data_array = zeros(ntime, ndoses, nextracts, nspecies);
data_array(:, :, 1, 1) = dG(:, CUSTOMCOLIDX{2,1});
data_array(:, :, 2, 1) = dG(:, CUSTOMCOLIDX{2,2});
data_array(:, :, 3, 1) = dG(:, CUSTOMCOLIDX{2,3});

% the order needed for the function log_likelihood_sharedCSP is 
%  nT x nMS x nIC x nEnv:
data_array = permute(data_array, [1, 4, 2, 3]); 

end
