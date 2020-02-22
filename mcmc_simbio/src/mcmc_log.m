function mcmc_log(tstamp, projdir,savedir, mi, di)
% mcmc_log generate mcmc log text file. 
%   This function generates a log file describing the mcmc run used to call
%   it. Inputs are the timestamp string: tstamp, the savedir string, ie, which 
% directory to save in, and the mcmc_info struct, mi. 
% The log file has information on:
% 1. Any descriptive info about the estimation run specified in the mcmc_info
% struct. 
% 2. The name of the model, if one is present in the model object. 
% 3. MCMC parameters: path, stdev, # of walkers, step size, tightening,
% number of repeats, thinnning, points per iteration, number of mcmc steps
% (~= points per iteration / (thinning * number of walkers)
% 4. List of estimated parameters, and the parameter ranges of the
% parameters
% 5. Experimental setup: dosed species, measured species (in additive groups),
% 
% FUTURE VERSIONS:
% 6. Model equations
% 7. Simulation status: Completed without errors, Incomplete due to user 
% termination, Incomplete for other reason. If incomplete, print out the 
% error stack
% 8.  A list of up to 50 sample points from the parmeter vector
% 
% !TODO (later) generalize this to the multi model mode. 
% 

currdir = pwd;
cd(savedir);
fileID = fopen(['summary_' tstamp '.txt'],'w');

%% Circuit description
if ~isempty(mi.circuitInfo)
    fS = ['\n ################################################ \n' ...
        'Circuit description: \n'];
    fprintf(fileID,fS);
    fS = mi.circuitInfo; 
    % circuit info should be a preformatted string for printing with fprintf
    fprintf(fileID,fS);
end

%% Data Description
if ~isempty(di.dataInfo)
fS = ['\n ################################################ \n' ...
    'Data description: \n'];
fprintf(fileID,fS);
% data info should be a preformatted string for printing with fprintf
fS = di.dataInfo;
fprintf(fileID,fS);
end

if ~isempty(mi.modelName)
%% Model Name
fS = ['\n ################################################ \n' ...
    'model name: \n'];
fprintf(fileID,fS);
fS = mi.modelName;
fprintf(fileID,fS);
end

%% Simulation Parameters
fS = ['################################################ \n'...
    'Simulation Parameters are: \n'];

fprintf(fileID,fS);
fprintf(fileID,'path: %s \n', projdir);
fprintf(fileID,'stdev: %3.1f \n', mi.stdev);
fprintf(fileID,'number of walkers: %d \n', mi.nW);
fprintf(fileID,'step size:%4.2f \n', mi.stepSize);
fprintf(fileID,'tightening: %4.2f \n', mi.tightening);
fprintf(fileID,'number of repeats: %d \n', mi.nIter);
fprintf(fileID,'thinning: %d \n', mi.thinning);
fprintf(fileID,'points per iter: %d \n', mi.nPoints);

%% Estimated Parameter Info
fS = ['\n ################################################ \n' ...
    'The estimated species and parameters are: \n'];
fprintf(fileID,fS);
fS = '%s in range [%0.5g %0.5g] \n';
[nn] = length(mi.names_ord);
for i = 1:nn
    fprintf(fileID,fS,mi.names_ord{i}, (mi.paramRanges(i, 1)),...
        (mi.paramRanges(i, 2)));
end

%% Dosed Species Info
% names of dosed species from the mcmc_info struct and the data_info struct. 
fS = ['################################################ \n' ...
    'Dosed Species are (name specified in mcmc_info, name specified in data_info) \n'];
fprintf(fileID,fS);

dN1 = mi.dosedNames;
dN2 = di.dosedNames;
nn1 = length(dN1);
nn2 = length(dN2);

if nn1 ~= nn2 
    warning('the dosed names in the mcmc_info and the data_info structs must be the same.')
    fS = '%s \n';
    for i = 1:nn1
        fprintf(fileID,fS,dN1{i});
    end
else
    fS = '%s      and      %s \n';
    for i = 1:nn1
        fprintf(fileID,fS,dN1{i}, dN2{i});
    end
end



% a matrix specifying the set of dose combinations used
fS = ['################################################ \n' ...
    'Dosed combinations (in #species by #combinations) are: \n'];
fprintf(fileID,fS);
doseMatrix = mi.dosedVals;

for ii = 1:size(doseMatrix,1)
    fS = '%s \n';
    fprintf(fileID,fS,dN1{ii});
    fprintf(fileID,'%g\t',doseMatrix(ii,:));
    fprintf(fileID,'\n');
end



fS = ['\n ################################################ '...
'\n The Measured Species are: \n'];
fprintf(fileID,fS);

fS = '%s \n';
fS1 = 'Measured Species group number %d comprises of: \n';
[nn] = length(mi.measuredSpecies); % number of species groups
for i = 1:nn
	fprintf(fileID,fS1, i);
	speciesList = mi.measuredSpecies{i};

	for j = 1:length(speciesList)
    fprintf(fileID,fS,speciesList{j});
	end

	fprintf(fileID, '--------------------------------------------- \n');

end

fclose(fileID);
cd(currdir)

end

