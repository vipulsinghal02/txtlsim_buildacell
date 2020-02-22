function mcmc_log_v2(tstamp,projdir,  savedir, mcmc_info, di, initialization_used, stepMultipliers, prevtstamp, temperature)
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
ri = mcmc_info.runsim_info;
mi = mcmc_info.model_info;
mai = mcmc_info.master_info;
currdir = pwd;
cd(savedir);
fileID = fopen(['summary_' tstamp '.txt'],'w');

initialization_used
%% Simulation Parameters
fS = ['################################################ \n'...
    'Simulation Parameters are: \n'];

fprintf(fileID,fS);
fprintf(fileID,'path: %s \n', projdir);
fprintf(fileID,'stdev: %3.1f \n', ri.stdev);
fprintf(fileID,'number of walkers: %d \n', ri.nW);
fprintf(fileID,'step size:%4.2f \n', ri.stepSize);
fprintf(fileID,'tightening: %4.2f \n', ri.tightening);
fprintf(fileID,'number of repeats: %d \n', ri.nIter);
fprintf(fileID,'thinning: %d \n', ri.thinning);
fprintf(fileID,'points per iter: %d \n', ri.nPoints);
fprintf(fileID,['step size ladder was: ' num2str(stepMultipliers) '\n']);
fprintf(fileID,'prevtstamp: %s \n', prevtstamp);

for i = 1:length(mi)
    %% Circuit description
    if  ~isempty(mi(i).circuitInfo)
        fS = ['\n ################################################ \n' ...
            'Circuit description: \n'];
        fprintf(fileID,fS);
           fS = [mi(i).circuitInfo '\n']; 
            % circuit info should be a preformatted string for printing with fprintf
            fprintf(fileID,fS); 
    end

    if  ~isempty(mi(i).modelName)
        %% Model Name
        fS = ['\n ################################################ \n' ...
            'model name: \n'];
        fS = [mi(i).modelName '\n'];
        fprintf(fileID,fS);
    end


    %% Estimated Parameter Info
    fS = ['\n ################################################ \n' ...
        'The estimated species and parameters are (log transformed): \n'];
    fprintf(fileID,fS);
    fS = '%s in range [%0.5g %0.5g] \n';
    [nn] = length(mai.estNames);
    for j = 1:nn
        fprintf(fileID,fS,mai.estNames{j}, (mai.paramRanges(j, 1)),...
            (mai.paramRanges(j, 2)));
    end

    fS = ['################################################ \n'...
        'initialization used was: \n'];
        fprintf(fileID,fS);
        fS = [initialization_used '\n'];
        fprintf(fileID,fS);
    % %% Dosed Species Info
    % % names of dosed species from the mcmc_info struct and the data_info struct. 
    % fS = ['################################################ \n' ...
    %     'Dosed Species are (name specified in mcmc_info, name specified in data_info) \n'];
    % fprintf(fileID,fS);

    fS = ['\n ################################################ '...
    '\n The Measured Species are: \n'];
    fprintf(fileID,fS);

    fS = '%s \n';
    fS1 = 'Measured Species group number %d comprises of: \n';
    [nn] = length(mi(i).measuredSpecies); % number of species groups
    for j = 1:nn
        fprintf(fileID,fS1, j);
        speciesList = mi(i).measuredSpecies{j};

        for k = 1:length(speciesList)
        fprintf(fileID,fS,speciesList{k});
        end

        fprintf(fileID, '--------------------------------------------- \n');

    end

end

%% Data Description
for i = 1:length(di)
    if  ~isempty(di(i).dataInfo)
        fS = ['\n ################################################ \n' ...
            'Data description: \n'];
        % data info should be a preformatted string for printing with fprintf
        fS = [di(i).dataInfo '\n'];
        fprintf(fileID,fS);
    end
    
end


fclose(fileID);
cd(currdir);





% dN1 = mi.dosedNames;
% dN2 = di.dosedNames;
% nn1 = length(dN1);
% nn2 = length(dN2);

% if nn1 ~= nn2 
%     warning('the dosed names in the mcmc_info and the data_info structs must be the same.')
%     fS = '%s \n';
%     for i = 1:nn1
%         fprintf(fileID,fS,dN1{i});
%     end
% else
%     fS = '%s      and      %s \n';
%     for i = 1:nn1
%         fprintf(fileID,fS,dN1{i}, dN2{i});
%     end
% end



% % a matrix specifying the set of dose combinations used
% fS = ['################################################ \n' ...
%     'Dosed combinations (in #species by #combinations) are: \n'];

% doseMatrix = mi.dosedVals;

% for ii = 1:size(doseMatrix,1)
%     fS = '%s \n';
%     fprintf(fileID,fS,dN1{ii});
%     fprintf(fileID,'%g\t',doseMatrix(ii,:));
%     fprintf(fileID,'\n');
% end





end

