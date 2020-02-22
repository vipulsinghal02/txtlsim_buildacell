function mcmc_runsim(tstamp, projdir, tv,da, em, mi)
%mcmc_runsim run the mcmc estimation 
% em is the exported simbiology model
% mi is the mcmc info struct created using the mcmc_info_... function. 
% funcs are the various functions needed, the prior, the lognormvec. 

% !LATER check if a directory to save the simulation results exists, and create it
% if it does not. 
% for now we just assume that project_init did its job right. 


% order the est names in the same order as 
% run an initial pass to find a set of points where integration tolerances
% are met. This valid points will be our initial walker positions. 
eno = mi.names_ord;
ds = struct('names', {mi.dosednames}, 'dosematrix', mi.dosedvals);
minit = integrableLHS(em, mi.nW, mi.paramranges, eno, ...
    ds);

% setup the log prior, log likelihood function and lognormvec functions
lognormvec=@(res,sig) -(res./sig).^2 -log(sqrt(2*pi)).*sig;

logprior = @(logp) all(mi.paramranges(:, 1) < logp) &&...
    all(logp < mi.paramranges(:,2));

sigg = mi.stdev/mi.tightening;

% in the future, change this to use input supplied functions in place of
% lognormvec and gen_residuals_3. 
loglike = @(logp) sum(lognormvec(gen_residuals_3(logp, em, da, tv, ...
    mi.dosedvals, mi.measuredspecies),sigg));

% run the burn in simulation
 
        tic
        [m] =gwmcmc_vse(minit,{logprior loglike},...
            mi.npoints,...
            'StepSize',mi.stepsize , ...
            'ThinChain',mi.thinning,...
            'Parallel', mi.parallel);
        toc
        
%         pause(10)
        
        minit = m(:,:,end);
        clear m
    % run the actual simuations, saving the data every iteration
    cfname = cell(mi.niter, 1);
    specificproj = [projdir '/simdata_' tstamp];
        for i = 1:mi.niter %
            tic
            disp(sprintf('starting mcmc %d\n', i));
             [m, ~, ~, ~] = gwmcmc_vse(minit,{logprior loglike},...
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

siminfo = {'MG apt and GFP', '''data_dsg2014'''};
fS = 'MCMC estimation for %s, using the TXTL modeling toolbox \n';
fprintf(fileID,fS,siminfo{1});
fS = 'Data from file %s \n';
fprintf(fileID,fS, siminfo{2});

fS = '################################################ \nThe estimated species and parameters are: \n';
fprintf(fileID,fS);
fS = '%s in range [%0.5g %0.5g] \n';
[nn] = length(eno);
for i = 1:nn
    fprintf(fileID,fS,eno{i}, exp(mi.paramranges(i, 1)),...
                  exp(mi.paramranges(i, 2)));
end

fS = '################################################ \nDosed Species are: \n';
fprintf(fileID,fS);
fS = '%s \n';
[nn] = length(ds.names);
for i = 1:nn
    fprintf(fileID,fS,ds.names{i});
end

% do the measured species later: 
% fS = '################################################ \nMeasured Species are: \n';
% fprintf(fileID,fS);
% fS = '%s \n';
% [nn] = length(mi.measuredspecies);
% for i = 1:nn
%     fprintf(fileID,fS,mi.measuredspecies{i});
% end

fS = '################################################ \nSimulation Parameters are: \n';
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

