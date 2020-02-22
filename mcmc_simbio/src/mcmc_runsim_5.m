function mcmc_runsim_5(tstamp, projdir, tv,da, em, mi, pmap, varargin)
%mcmc_runsim run the mcmc estimation - This version works with (the
%prototyping version) the gen_residuals_5 function. 
% 
% em is the exported simbiology model
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
% order the est names in the same order as 
% run an initial pass to find a set of points where integration tolerances
% are met. This valid points will be our initial walker positions. 

p = inputParser;
p.addParameter('distribution', 'LHS', @ischar)
p.addParameter('width', 0.1, @isnumeric)
p.addParameter('userinitialize', [], @isnumeric)

p.parse(varargin{:});

p = p.Results;
eno = mi.names_ord;
ds = struct('names', {mi.dosednames}, 'dosematrix', mi.dosedvals);

if isempty(p.userinitialize)
    % need to generalize this to the multi geometry case. 
    minit = integrableLHS(em, mi.nW, mi.paramranges, eno, ...
        ds, 'distribution', p.distribution, 'width', p.width);
else
    minit = p.userinitialize; 
    % assume all the user defined points are integrable. 
end

% setup the log prior, log likelihood function and lognormvec functions
lognormvec=@(res,sig) -(res./sig).^2 -log(sqrt(2*pi)).*sig;

espix = pmap{1};
cspix = pmap{2};
nExt = size(da, 5);


% nparam = nExt*length(espix) + length(cspix);
paramranges = mi.paramranges;

priorboundz = [repmat(paramranges(espix,:), length(espix),1) ; paramranges(cspix,:)];


logprior = @(logp) all(priorboundz(:, 1) < logp) &&...
    all(logp < priorboundz(:,2));

sigg = mi.stdev/mi.tightening;

loglike = @(logp) gen_residuals_5(logp, em, da, tv, ...
    mi.dosedvals, mi.measuredspecies, lognormvec, sigg, pmap );


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

    
    % run the actual simuations, saving the data every iteration
    cfname = cell(mi.niter, 1);
    specificproj = [projdir '/simdata_' tstamp];
        for i = 1:mi.niter %
            pause(10)
            tic
            disp(sprintf('starting mcmc %d\n', i));
             [m] = gwmcmc_vse(minit,{logprior loglike},...
                 mi.npoints, ...
                 'StepSize',mi.stepsize , ...
            'ThinChain',mi.thinning, 'Parallel', mi.parallel);
        
            disp(sprintf('ending mcmc %d\n', i));
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


%% we save all the other things here (ie, outside the loops
fname = ['full_variable_set_' tstamp]; % filename
save([specificproj '/' fname]);



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


% !TODO : THIS IS NOT RIGHT. 
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

