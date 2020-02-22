% Vipul Singhal
% Caltech Dec 2014 (edits: 2019)

%% Import data and convert to groupedData object for fitting
 
txtlDir = '/Users/vipulsinghal/Dropbox/Documents/toolbox/txtlsim_vsfork2017';

cd(txtlDir)
addpath(txtlDir);
txtl_init
mcmc_init


addpath(genpath('/Users/vipulsinghal/Dropbox/110114/Vipul_computational/Computational Tools'));
datadir = '/Users/vipulsinghal/Dropbox/110114/Vipul_computational/Computational Tools/RawComputationalData';
load([datadir '/Set6Data.mat'],'Set6Data')
load([datadir '/Set7Data.mat'],'Set7Data')
load([datadir '/Set8Data.mat'],'Set8Data')
load([datadir '/Set9Data.mat'],'Set9Data')
load([datadir '/Set10Data.mat'],'Set10Data')

ID = reshape(repmat(linspace(1,35,35), 81, 1), 81*35,1);

time = [repmat(Set6Data(:,1,1),7,1);
    repmat(Set7Data(:,1,1),7,1)
    repmat(Set8Data(:,1,1),7,1);
    repmat(Set9Data(:,1,1),7,1)
    repmat(Set10Data(:,1,1),7,1);];

GFP  = 1000*[mean(Set6Data(:,2:4,1),2); 
    mean(Set6Data(:,2:4,2),2);
    mean(Set6Data(:,2:4,3),2);
    mean(Set6Data(:,2:4,4),2);
    mean(Set6Data(:,2:4,5),2);
    mean(Set6Data(:,2:4,6),2);
    mean(Set6Data(:,2:4,7),2);
    mean(Set7Data(:,2:4,1),2); 
    mean(Set7Data(:,2:4,2),2);
    mean(Set7Data(:,2:4,3),2);
    mean(Set7Data(:,2:4,4),2);
    mean(Set7Data(:,2:4,5),2);
    mean(Set7Data(:,2:4,6),2);
    mean(Set7Data(:,2:4,7),2);
    mean(Set8Data(:,2:3,1),2); 
    mean(Set8Data(:,2:4,2),2);
    mean(Set8Data(:,2:4,3),2);
    mean(Set8Data(:,2:4,4),2);
    mean(Set8Data(:,2:4,5),2);
    mean(Set8Data(:,2:4,6),2);
    mean(Set8Data(:,2:4,7),2);
    mean(Set9Data(:,2:4,1),2); 
    mean(Set9Data(:,2:4,2),2);
    mean(Set9Data(:,2:4,3),2);
    mean(Set9Data(:,2:4,4),2);
    mean(Set9Data(:,2:4,5),2);
    mean(Set9Data(:,2:4,6),2);
    mean(Set9Data(:,2:4,7),2);
    mean(Set10Data(:,2:4,1),2); 
    mean(Set10Data(:,2:4,2),2);
    mean(Set10Data(:,2:4,3),2);
    mean(Set10Data(:,2:4,4),2);
    mean(Set10Data(:,2:4,5),2);
    mean(Set10Data(:,2:4,6),2);
    mean(Set10Data(:,2:4,7),2);];

GFP = max(GFP,1e-10);

% tetR DNA
d3OC = NaN(35*81,1); %1000 nM default
dLASR = NaN(35*81,1); % 1nM
dATC = NaN(35*81,1); % 10000nM
dTETR = NaN(35*81,1); % 0.1nM
dGFP = NaN(35*81,1); % 1nM


d3OC(1:81:81*35) = 1000;
dLASR(1:81:81*35) = 1; % 1nM
dATC(1:81:81*35) = 10000; % 10000nM
dTETR(1:81:81*35) = 0.1; % 0.1nM
dGFP(1:81:81*35) = 1; % 1nM



conc3OC = [10000 1000 100 10 1 0.1 0]'; 
conclasR = [2 1 0.5 0.25 0.125 0.0625 0]';
concaTc = [10000 1000 100 10 1 0.1 0]';
conctetR = [1 0.1 0.01 0.001 0.0001 0.00001 0]';
concGFP = [4 2 1 0.5 0.25 0.125 0]';

for i = 1:7
    d3OC(81*(i-1)+1) = conc3OC(i);
    dLASR(81*(7+i-1)+1) = conclasR(i);
    dATC(81*(14+i-1)+1) = concaTc(i);
    dTETR(81*(21+i-1)+1) = conctetR(i);
    dGFP(81*(28+i-1)+1) = concGFP(i);
end

% set aTc to zero for the tetR DNA variation case (ID 22 - 28)
dATC(1+81*21:81:1+81*27) = 0;


DataTable = table(ID,time, GFP, d3OC, dLASR, dATC, dTETR, dGFP);

grpData_test = groupedData(DataTable);
% save('/Users/vipulsinghal/Dropbox/110114/Vipul_computational/Computational Tools/s6tos10_combined_grpData_corrected.mat', 'grpData_test');


% 
% m = Mobj_s3s4;
% load('/Users/vipulsinghal/Dropbox/110114/Vipul_computational/Computational Tools/EstimationTests/s3s4_pre_estParams_2015January06_104501.mat') 
% updateParam = {'kf',   0.2588;
% 'kr',327.298;
% 'ptet_sequestration_F', 77.67
% 'ptet_sequestration_R', 24.109
% 'TXTL_INDUCER_TETR_ATC_F', 0.397084;
% 'TXTL_INDUCER_TETR_ATC_R', 0.7782};
% 
% % Use this data to update the 'params' array:
% for i = 1:size(updateParam,1)
%     idx = find(strcmp(updateParam{i,1}, params(:,2)));
%     if  ~isempty(idx)
%         params{idx, 3} = updateParam{i,2};
%         
%     else
%         error('parameter not found in ''params'' array');
%     end
% end
% 
% dose = {'DNA placNov4--rbs--tetRNov4', [2 0.2 0.02 0.002 0.0002 0.00002 0 0.1 0.1 0.1 0.1 0.1 0.1 0.1];
%     'aTc', [0 0 0 0 0 0 0 10000 1000 100 10 1 0.1 0]};
% titl = 'GFP levels, plac-tetR varying, or aTc varying';
% legs = {'2nM', '0.2nM', '0.02nM', '0.002nM', '0.0002nM', '0.00002nM', '0nM', '10000nM', '1000nM', '100nM', '10nM', '1nM', '0.1nM'};
% extraInfo = {'1mM IPTG';'ptet-UTR1-deGFP 1nM';'plac-tetR varying'};
% saveDir = '/Users/vipulsinghal/Dropbox/110114/Vipul_computational/Computational Tools/EstimationTests/component fits/s3s4 iterations/';
% saveFile = ['s3s4_combined' datestr(now,'yyyymmmmdd_HHMMSS')];
% % txtlDir, m, estParam, dose, dataArray, titl, legs, extraInfo, saveDir, savefile
% % f = GenerateTrainingPlots(txtlDir, m, params, dose, Set3Data, titl, legs, extraInfo, saveDir, saveFile);
% 
% 
