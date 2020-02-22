%% Step 1: Build time vectors and data matrices
addpath(genpath('/Users/vipulsinghal/Dropbox/Documents/MATLAB'))
tv =  60*(0:8:(8*60))';
dmexp = zeros(length(tv), 3*35);% data matrix of experimental data in triplicate

% load S1 data
load('/Users/vipulsinghal/Dropbox/110114/Vipul_computational/Computational Tools/S1/S1data.mat');
iGFP = findspecies(m, 'protein deGFP*');
for i = 1:7
    % trajectories of original data in triplicate
    [T, it_exp, ~] = unique(Set1Data(:,1,i));
    dmexp(:,3*(i-1)+(1:3))= interp1(T, Set1Data(it_exp,2:4,i),tv);
end

% load S2 data 
load('/Users/vipulsinghal/Dropbox/110114/Vipul_computational/Computational Tools/S2/S2data.mat');
iGFP = findspecies(m, 'protein deGFP*');
for i = 1:7
    offset = 7;
    [T, it_exp, ~] = unique(Set2Data(:,1,i));
    dmexp(:,3*(offset+(i-1))+(1:3))= interp1(T, Set2Data(it_exp,2:4,i),tv);
end

% load S3 and S4 data 
load('/Users/vipulsinghal/Dropbox/110114/Vipul_computational/Computational Tools/S3/S3S4data.mat');
iGFP = findspecies(m, 'protein deGFP*');

t_ori = T_ORI{1}; % first element is S3
x_ori = X_ORI{1};
t_est = T_EST{1};
x_est = X_EST{1}; 

for i = 1:7
    offset = 14;
    [T, it_exp, ~] = unique(Set3Data(:,1,i));
    dmexp(:,3*(offset+(i-1))+(1:3))= interp1(T, Set3Data(it_exp,2:4,i),tv);
    
end

t_ori = T_ORI{2}; % second element is S4
x_ori = X_ORI{2};
t_est = T_EST{2};
x_est = X_EST{2}; 

iGFP = findspecies(m, 'protein deGFP*');
for i = 1:7
    offset = 21;
    [T, it_exp, ~] = unique(transformedSet4data(:,1,i));
    dmexp(:,3*(offset+(i-1))+(1:3))= interp1(T, transformedSet4data(it_exp,2:4,i),tv);
end

% load S5 data 
load('/Users/vipulsinghal/Dropbox/110114/Vipul_computational/Computational Tools/S5/S5data.mat');
iGFP = findspecies(m, 'protein deGFP*');
for i = 1:7
    offset = 28;
    % trajectories of original data in triplicate
    [T, it_exp, ~] = unique(Set5Data(:,1,i));
    dmexp(:,3*(offset+(i-1))+(1:3))= interp1(T, Set5Data(it_exp,2:4,i),tv);
end

%% Step 2: build index arrays

% experimental data
ixexp = {[1 2 3; 4 5 6; 7 8 9; 10 11 12; 13 14 15; 16 17 18; 19 20 21], [];
    14*3+[1 2 3; 4 5 6; 7 8 9; 10 11 12; 13 14 15; 16 17 18; 19 20 21], [];
    21*3+[1 2 3; 4 5 6; 7 8 9; 10 11 12; 13 14 15; 16 17 18; 19 20 21], [];
     7*3+[1 2 3; 4 5 6; 7 8 9; 10 11 12; 13 14 15; 16 17 18; 19 20 21], [];
    28*3+[1 2 3; 4 5 6; 7 8 9; 10 11 12; 13 14 15; 16 17 18; 19 20 21], []};

%% Step 3: build title arrays
textarray = {'pTet constitutive data', 'pTet constitutive fit'; 
    'TetR repression data','TetR repression fit'; 
    'aTc induction data','aTc induction fit';
    'pLac constitutive data','pLac constitutive fit';
    '3oc12 AHL induction data','3oc12 AHL induction fit'};
ta = struct('text', {textarray},...
    'fontsize', 16);

%% Step 4: Build legend arrays
textarray = {{'4nM', '2nM', '1nM', '0.5nM', '0.25nM', '0.125nM', '0nM'};
    {'2nM', '0.2nM', '0.02nM', '0.002nM', '0.0002nM', '0.00002nM', '0 nM'};
    {'10000nM', '1000nM', '100nM', '10nM', '1nM', '0.1nM', '0 nM'};
    {'4nM', '2nM', '1nM', '0.5nM', '0.25nM', '0.125nM', '0nM'};
    {'10000nM', '1000nM', '100nM', '10nM', '1nM', '0.1nM', '0 nM'}};
la = struct('text', {textarray},...
    'fontsize', 9, ...
    'location', 'NorthWest');


%% Step 5: Build ylabels
textarray = {'GFP, uM','GFP, uM'; 
    'GFP, uM','GFP, uM'; 
    'GFP, uM','GFP, uM'; 
    'GFP, uM','GFP, uM'; 
    'GFP, uM','GFP, uM'};
yl = struct('text', {textarray},...
    'fontsize', 12,...
    'axislimit', 'rowmax');% 'individual', 'rowmax', or a scalar input of type double


%% Step 5: Build xlabels
textarray = {'time, hours','time, hours'; 
    'time, hours','time, hours'; 
    'time, hours','time, hours'; 
    'time, hours','time, hours'; 
    'time, hours','time, hours'};
xl = struct('text', {textarray},...
    'fontsize', 12, ...
    'axislimit', 8);

%% Step 6: Build the data structures

scrsz = get(groot,'ScreenSize')
fs= [100 scrsz(4)/8 scrsz(3)/1.5 scrsz(4)/1.05]; %[left bottom width height]
% construct the relevant fields for the data structs
close all

expdata = struct('d', dmexp,...
    't', tv/3600,...
    'Ix', {ixexp},...
    'estimator', 'mean',...
    'errb', true,...
    'color', @parula,...
    'linespec', '--',...
    'linewidth',1,...
    'titles', {[]},...
    'legends', {[]},...
    'ylabel', {[]},...
    'xlabel',{[]},...
    'figuresize', {[]},...
    'userdata', {[]});

%% Step 7: Plot away
plot_multiple(expdata);

