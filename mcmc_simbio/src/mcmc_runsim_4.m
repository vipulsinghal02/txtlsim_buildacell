function mcmc_runsim_4(tstamp, projdir, tv,da, mobj, mi, varargin)
% mcmc_runsim run the mcmc estimation - This version works with (the
% prototyping version) the gen_residuals_4 function.
%
% mobj is the simbiology model
% mi is the mcmc info struct created using the mcmc_info_... function.
% funcs are the various functions needed, the prior, the lognormvec.

% !LATER check if a directory to save the simulation results exists,
% and create it if it does not.
% for now we just assume that project_init did its job right.
%
% has the following name value input pairs.
% 'distribution', 'LHS'
% 'width', 0.1
% 
% 
% 
% order the est names in the same order as
% run an initial pass to find a set of points where integration tolerances
% are met. This valid points will be our initial walker positions.

p = inputParser;
p.addParameter('distribution', 'LHS', @ischar)
p.addParameter('width', 0.1, @isnumeric)
p.addParameter('userinitialize', [], @isnumeric)

p.parse(varargin{:});

p = p.Results;

%% EXPORT MODEL object to get it ready for MCMC
% the resulting object is of class SimBiology.export.Model
% documentation: https://www.mathworks.com/help/simbio/ref/simbiology.export.model-class.html
% Sven Mesecke's blog post on using the exported model class for parameter inference applicaton. 
% http://sveme.org/how-to-use-global-optimization-toolbox-algorithms-for-simbiology-parameter-estimation-in-parallel-part-i.html

% export and accelerate simbiology model object using estimated species
% and dosing species names

% select the parameters and species objects using the name array
ep = sbioselect(mobj, 'Type', 'parameter', 'Name', ...
    mi.names_unord);% est parameters

es = sbioselect(mobj, 'Type', 'species', 'Name', ...
    mi.names_unord);% est species

aps = [ep; es]; % active parameters and species

% reorder the parameter and species so they are in the same order as that
% in the model. 
eno = cell(length(aps), 1);% est names ordered

for i = 1:length(aps)
    eno{i} = aps(i).Name;
end

% 
ds = sbioselect(mobj, 'Type', 'species', 'Name', mi.dosednames);

emo = export(mobj, [ep; es; ds]); % exported model object, dosed species names. 
SI = emo.SimulationOptions;
SI.StopTime = tv(end);
accelerate(emo);

mi.names_ord = eno; % est names ordered. 
mi.emo = emo; % exported model object. 


ds = struct('names', {mi.dosednames}, 'dosematrix', mi.dosedvals);
if isempty(p.userinitialize)
    minit = integrableLHS(em, mi.nW, mi.paramranges, eno, ...
        ds, 'distribution', p.distribution, 'width', p.width);
else
    minit = p.userinitialize;
    % assume all the user defined points are integrable.
end

% setup the log prior, log likelihood function and lognormvec functions
lognormvec=@(res,sig) -(res./sig).^2 -log(sqrt(2*pi)).*sig;

logprior = @(logp) all(mi.paramranges(:, 1) < logp) &&...
    all(logp < mi.paramranges(:,2));

sigg = mi.stdev/mi.tightening;

loglike = @(logp) gen_residuals_4(logp, em, da, tv, ...
    mi.dosedvals, mi.measuredspecies, lognormvec, sigg);


% run the burn in simulation
if isempty(p.userinitialize)
    tic
    [m] =gwmcmc_vse(minit,{logprior loglike},...
        mi.npoints,...
        'StepSize',mi.stepsize , ...
        'ThinChain',mi.thinning,...
        'Parallel', mi.parallel);
    toc
    
    minit = m(:,:,end);
    clear m
else
    disp('User initialized intitial walker positions, skipping burn in phase')
end


% specify where to save things 
cfname = cell(mi.niter, 1);
specificproj = [projdir '/simdata_' tstamp];

%% we save useful variables in a one off manner here (ie, outside the loops)
fname = ['full_variable_set_' tstamp]; % filename
save([specificproj '/' fname]);

% run the actual simuations, saving the data every iteration
for i = 1:mi.niter %
    if ~mod(i, 3)
        fprintf('Pausing for 10 minutes before starting run number %d. \n', i);
        pause(600)        
    end
    
    tic
    fprintf('starting mcmc %d\n', i);
    [m] = gwmcmc_vse(minit,{logprior loglike},...
        mi.npoints, ...
        'StepSize',mi.stepsize , ...
        'ThinChain',mi.thinning, 'Parallel', mi.parallel);
    
    fprintf('ending mcmc %d\n', i);
    toc
    fname = ['mcmc' tstamp '_ID' num2str(i)] ;
    cfname{i} = fname;
    save([specificproj '/' fname], 'm');
    % the only thing that is different in each run are the above
    %             pause(1)
    minit = m(:,:,end);% + 0.1*randn(size(m(:,:,end-1)));
    
    clear m
    
end

% generate log file



%% generate simulation log file
% do this. also, write a log file using log4m and
%  fprintf.
currdir = pwd;
cd(specificproj);
fileID = fopen(['summary_' tstamp '.txt'],'w');

% !TODO FIX THIS. REALLY SHOULD NOT BE LIKE THIS.
% siminfo = {'MG apt and GFP', '''data_dsg2014'''};
% fS = 'MCMC estimation for %s, using the TXTL modeling toolbox \n';
% fprintf(fileID,fS,siminfo{1});
% fS = 'Data from file %s \n';
% fprintf(fileID,fS, siminfo{2});

fS = ['################################################ \n' ...
    'The estimated species and parameters are: \n'];
fprintf(fileID,fS);
fS = '%s in range [%0.5g %0.5g] \n';
[nn] = length(eno);
for i = 1:nn
    fprintf(fileID,fS,eno{i}, exp(mi.paramranges(i, 1)),...
        exp(mi.paramranges(i, 2)));
end

fS = ['################################################ \n' ...
    'Dosed Species are: \n'];
fprintf(fileID,fS);
fS = '%s \n';
[nn] = length(ds.names);
for i = 1:nn
    fprintf(fileID,fS,ds.names{i});
end

% do the measured species later:
% fS = ['################################################ '...
% '\nMeasured Species are: \n'];
% fprintf(fileID,fS);
% fS = '%s \n';
% [nn] = length(mi.measuredspecies);
% for i = 1:nn
%     fprintf(fileID,fS,mi.measuredspecies{i});
% end

fS = ['################################################ \n'...
    'Simulation Parameters are: \n'];
fprintf(fileID,fS);
fprintf(fileID,'path: %s \n', projdir);
fprintf(fileID,'stdev: %3.1f \n', mi.stdev);
fprintf(fileID,'number of walkers: %d \n', mi.nW);
fprintf(fileID,'step size:%4.2f \n', mi.stepsize);
fprintf(fileID,'tightening: %4.2f \n', mi.tightening);
fprintf(fileID,'number of repeats: %d \n', mi.niter);
fprintf(fileID,'thinning: %d \n', mi.thinning);
fprintf(fileID,'points per iter: %d \n', mi.npoints);

% MAP estimates
% median
fclose(fileID);
cd(currdir)

end

