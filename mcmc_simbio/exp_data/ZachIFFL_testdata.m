function di = ZachIFFL_testdata(varargin)

% Vipul Singhal
% Caltech Dec 2014 (edits: 2019)

%% Import data and convert to groupedData object for fitting

datamode = 'means';
if nargin == 1
    datamode = varargin{1};
    if ~strcmp(datamode, 'means') && ~strcmp(datamode, 'all_trajectories')
        error('Not a valid data mode.')
    end
end

txtlDir = '/Users/vipulsinghal/Dropbox/Documents/toolbox/txtlsim_vsfork2017';

cd(txtlDir)
addpath(txtlDir);
txtl_init
mcmc_init

datadir = '/Users/vipulsinghal/Dropbox/Documents/toolbox/txtlsim_vsfork2017/mcmc_simbio/exp_data/Zach_IFFL_raw/S6';
load([datadir '/Set6Data.mat'],'Set6Data')
load([datadir '/Set7Data.mat'],'Set7Data')
load([datadir '/Set8Data.mat'],'Set8Data')
load([datadir '/Set9Data.mat'],'Set9Data')
load([datadir '/Set10Data.mat'],'Set10Data')

if strcmp(datamode, 'means')
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
    % /Users/vipulsinghal/Dropbox/Documents/toolbox/txtlsim_vsfork2017/mcmc_simbio/exp_data/Zach_IFFL_raw/S6/s6tos10_combined_grpData_corrected.mat
    % grpData_test(1:81:35*81,:)
    %
    % ans =
    %
    %   35×8 table
    %
    %     ID    time      GFP      d3OC     dLASR     dATC     dTETR     dGFP
    %     __    ____    _______    _____    ______    _____    ______    _____
    %
    %      1     0        1e-10    10000         1    10000       0.1        1
    %      2     0        1e-10     1000         1    10000       0.1        1
    %      3     0       1.7691      100         1    10000       0.1        1
    %      4     0       2.6537       10         1    10000       0.1        1
    %      5     0        1e-10        1         1    10000       0.1        1
    %      6     0       2.5063      0.1         1    10000       0.1        1
    %      7     0      0.29485        0         1    10000       0.1        1
    %      8     0        1e-10     1000         2    10000       0.1        1
    %      9     0      0.29485     1000         1    10000       0.1        1
    %     10     0      0.29485     1000       0.5    10000       0.1        1
    %     11     0      0.14743     1000      0.25    10000       0.1        1
    %     12     0        1e-10     1000     0.125    10000       0.1        1
    %     13     0       4.4228     1000    0.0625    10000       0.1        1
    %     14     0        1.032     1000         0    10000       0.1        1
    %     15     0        1e-10     1000         1    10000       0.1        1
    %     16     0       1.7691     1000         1     1000       0.1        1
    %     17     0       3.8331     1000         1      100       0.1        1
    %     18     0       3.5383     1000         1       10       0.1        1
    %     19     0       5.3074     1000         1        1       0.1        1
    %     20     0       5.7497     1000         1      0.1       0.1        1
    %     21     0       5.8971     1000         1        0       0.1        1
    %     22     0       3.9805     1000         1        0         1        1
    %     23     0        4.128     1000         1        0       0.1        1
    %     24     0       4.5702     1000         1        0      0.01        1
    %     25     0       5.7497     1000         1        0     0.001        1
    %     26     0       3.2434     1000         1        0    0.0001        1
    %     27     0       4.4228     1000         1        0     1e-05        1
    %     28     0       1.4743     1000         1        0         0        1
    %     29     0        3.096     1000         1    10000       0.1        4
    %     30     0       3.5383     1000         1    10000       0.1        2
    %     31     0       4.8651     1000         1    10000       0.1        1
    %     32     0       4.2754     1000         1    10000       0.1      0.5
    %     33     0         5.16     1000         1    10000       0.1     0.25
    %     34     0       3.5383     1000         1    10000       0.1    0.125
    %     35     0        4.128     1000         1    10000       0.1        0
    
    %% build info that will go into the data info structs
    
    tv =  60*(0:8:(8*60))';
    dmexp = zeros(length(tv), 3*35);% data matrix of experimental data in triplicate
    
    datadescription = {...
        'Zach IFFL test data: vary 3OC12';
        'Zach IFFL test data: vary lasR DNA';
        'Zach IFFL test data: vary aTc';
        'Zach IFFL test data: vary tetR DNA';
        'Zach IFFL test data: vary deGFP DNA'};
    % final order is time x species x replicates x doses
    % da = permute(da, [1, 3, 4, 2]);
    dimlabels = {'time points', 'measured species', 'replicates', 'doses'};
    
    dv1 = [10000    1000    100     10      1       0.1     0];
    dv2 = [2        1       .5      .25     .125    .0625   .03125];
    dv3 = [10000    1000    100     10      1       0.1     0];
    dv4 = [1        0.1    0.01    0.001   0.0001   0.00001 0];
    dv5 = [4        2        1       .5      .25     .125   0];
    dv = {dv1, dv2, dv3, dv4, dv5};
    mn = {'deGFP'};
    dosedNames1 = {'3OC12HSL'};
    dosedNames2 = {'pLac-lasR DNA'};
    dosedNames3 = {'aTc'};
    dosedNames4 = {'pLas-tetR DNA'};
    dosedNames5 = {'pLas-tetO-deGFP'};
    dosedNames = {dosedNames1, dosedNames2, dosedNames3, dosedNames4, dosedNames5};
    
    tv = 0:360:28800 ;
    %% build the data array from the grouped data array.
    da = cell(5, 1);
    for i = 1:5
        da{i} = zeros(length(tv), 1, 1, 7);
        %     for rep = 1:3
        for doseid = 1:7
            currentIndexes = ((i-1)*(81*7)+((doseid-1)*81)+(1:81));
            da{i}(:, 1, 1, doseid) = grpData_test.GFP(currentIndexes);
        end
        %     end
    end
    
    %%
    % di = cell(5, 1)
    clear di
    for i = 1:5
        di(i) = struct('dataInfo', datadescription(i), ...
            'timeVector', {tv}, ...
            'timeUnits', {'seconds'},...
            'dataArray', da(i),...
            'measuredNames', {mn},...
            'dataUnits', {{'nM'}},...
            'dimensionLabels', {dimlabels}, ...
            'dosedNames', dosedNames(i),...
            'dosedVals', dv(i), ...
            'doseUnits', 'nM');
    end
    %%
    figure;
    for i = 1:5
        subplot(5, 1, i);
        plot(0:360/3600:28800/3600, squeeze(di(i).dataArray));
        title(di(i).dosedNames)
        
    end
    
elseif strcmp(datamode, 'all_trajectories')
    % somehow the set 8 data only has two columns. the third column is all
    % 0s . Dont remember why, but it is what it is.
    
%     Sets 6 and 7 have 3 replicates. 
% set 8 dose 1 has 2 replicates
% set 9 and 10 have 3 replicates. 
% Actually it is just easier to create a replicate. I know its not exactly
% right (it modifies the standard deviation), but it works with the data
% structures, adn does not affect the results in any substantial way. 

    
    
     
    GFP  = 1000*[Set6Data(:,2:4,1);
        Set6Data(:,2:4,2);
        Set6Data(:,2:4,3);
        Set6Data(:,2:4,4);
        Set6Data(:,2:4,5);
        Set6Data(:,2:4,6);
        Set6Data(:,2:4,7);
        Set7Data(:,2:4,1);
        Set7Data(:,2:4,2);
        Set7Data(:,2:4,3);
        Set7Data(:,2:4,4);
        Set7Data(:,2:4,5);
        Set7Data(:,2:4,6);
        Set7Data(:,2:4,7);
        [Set8Data(:,2:3,1) mean(Set8Data(:,2:3,1), 2)];
        Set8Data(:,2:4,2);
        Set8Data(:,2:4,3);
        Set8Data(:,2:4,4);
        Set8Data(:,2:4,5);
        Set8Data(:,2:4,6);
        Set8Data(:,2:4,7);
        Set9Data(:,2:4,1);
        Set9Data(:,2:4,2);
        Set9Data(:,2:4,3);
        Set9Data(:,2:4,4);
        Set9Data(:,2:4,5);
        Set9Data(:,2:4,6);
        Set9Data(:,2:4,7);
        Set10Data(:,2:4,1);
        Set10Data(:,2:4,2);
        Set10Data(:,2:4,3);
        Set10Data(:,2:4,4);
        Set10Data(:,2:4,5);
        Set10Data(:,2:4,6);
        Set10Data(:,2:4,7)];
    
    
    
    
    ID = reshape(repmat(linspace(1,35,35), 81, 1), 81*35,1);
    
    time = [repmat(Set6Data(:,1,1),7,1);
        repmat(Set7Data(:,1,1),7,1)
        repmat(Set8Data(:,1,1),7,1);
        repmat(Set9Data(:,1,1),7,1)
        repmat(Set10Data(:,1,1),7,1);];
    

    
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
    % /Users/vipulsinghal/Dropbox/Documents/toolbox/txtlsim_vsfork2017/mcmc_simbio/exp_data/Zach_IFFL_raw/S6/s6tos10_combined_grpData_corrected.mat
    % grpData_test(1:81:35*81,:)
    %
    % ans =
    %
    %   35×8 table
    %
    %     ID    time      GFP      d3OC     dLASR     dATC     dTETR     dGFP
    %     __    ____    _______    _____    ______    _____    ______    _____
    %
    %      1     0        1e-10    10000         1    10000       0.1        1
    %      2     0        1e-10     1000         1    10000       0.1        1
    %      3     0       1.7691      100         1    10000       0.1        1
    %      4     0       2.6537       10         1    10000       0.1        1
    %      5     0        1e-10        1         1    10000       0.1        1
    %      6     0       2.5063      0.1         1    10000       0.1        1
    %      7     0      0.29485        0         1    10000       0.1        1
    %      8     0        1e-10     1000         2    10000       0.1        1
    %      9     0      0.29485     1000         1    10000       0.1        1
    %     10     0      0.29485     1000       0.5    10000       0.1        1
    %     11     0      0.14743     1000      0.25    10000       0.1        1
    %     12     0        1e-10     1000     0.125    10000       0.1        1
    %     13     0       4.4228     1000    0.0625    10000       0.1        1
    %     14     0        1.032     1000         0    10000       0.1        1
    %     15     0        1e-10     1000         1    10000       0.1        1
    %     16     0       1.7691     1000         1     1000       0.1        1
    %     17     0       3.8331     1000         1      100       0.1        1
    %     18     0       3.5383     1000         1       10       0.1        1
    %     19     0       5.3074     1000         1        1       0.1        1
    %     20     0       5.7497     1000         1      0.1       0.1        1
    %     21     0       5.8971     1000         1        0       0.1        1
    %     22     0       3.9805     1000         1        0         1        1
    %     23     0        4.128     1000         1        0       0.1        1
    %     24     0       4.5702     1000         1        0      0.01        1
    %     25     0       5.7497     1000         1        0     0.001        1
    %     26     0       3.2434     1000         1        0    0.0001        1
    %     27     0       4.4228     1000         1        0     1e-05        1
    %     28     0       1.4743     1000         1        0         0        1
    %     29     0        3.096     1000         1    10000       0.1        4
    %     30     0       3.5383     1000         1    10000       0.1        2
    %     31     0       4.8651     1000         1    10000       0.1        1
    %     32     0       4.2754     1000         1    10000       0.1      0.5
    %     33     0         5.16     1000         1    10000       0.1     0.25
    %     34     0       3.5383     1000         1    10000       0.1    0.125
    %     35     0        4.128     1000         1    10000       0.1        0
    
    %% build info that will go into the data info structs
    
    tv =  60*(0:8:(8*60))';
    dmexp = zeros(length(tv), 3*35);% data matrix of experimental data in triplicate
    
    datadescription = {...
        'Zach IFFL test data: vary 3OC12';
        'Zach IFFL test data: vary lasR DNA';
        'Zach IFFL test data: vary aTc';
        'Zach IFFL test data: vary tetR DNA';
        'Zach IFFL test data: vary deGFP DNA'};
    % final order is time x species x replicates x doses
    % da = permute(da, [1, 3, 4, 2]);
    dimlabels = {'time points', 'measured species', 'replicates', 'doses'};
    
    dv1 = [10000    1000    100     10      1       0.1     0];
    dv2 = [2        1       .5      .25     .125    .0625   .03125];
    dv3 = [10000    1000    100     10      1       0.1     0];
    dv4 = [1        0.1    0.01    0.001   0.0001   0.00001 0];
    dv5 = [4        2        1       .5      .25     .125   0];
    dv = {dv1, dv2, dv3, dv4, dv5};
    mn = {'deGFP'};
    dosedNames1 = {'3OC12HSL'};
    dosedNames2 = {'pLac-lasR DNA'};
    dosedNames3 = {'aTc'};
    dosedNames4 = {'pLas-tetR DNA'};
    dosedNames5 = {'pLas-tetO-deGFP'};
    dosedNames = {dosedNames1, dosedNames2, dosedNames3, dosedNames4, dosedNames5};
    
    tv = 0:360:28800 ;
    %% build the data array from the grouped data array.
    da = cell(5, 1);
    for i = 1:5
        da{i} = zeros(length(tv), 1, 3, 7);
        %     for rep = 1:3
        for doseid = 1:7
            for repid = 1:3
                currentIndexes = ((i-1)*(81*7)+((doseid-1)*81)+(1:81));
                da{i}(:, 1, repid, doseid) = grpData_test.GFP(currentIndexes,repid);
            end
            
            
            
        end
        %     end
    end
    
    %%
    % di = cell(5, 1)
    clear di
    for i = 1:5
        di(i) = struct('dataInfo', datadescription(i), ...
            'timeVector', {tv}, ...
            'timeUnits', {'seconds'},...
            'dataArray', da(i),...
            'measuredNames', {mn},...
            'dataUnits', {{'nM'}},...
            'dimensionLabels', {dimlabels}, ...
            'dosedNames', dosedNames(i),...
            'dosedVals', dv(i), ...
            'doseUnits', 'nM');
    end
    
    
end

end

