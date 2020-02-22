function  mcmc_init
    % Initializes required paths. Call at the start of every new MATLAB session, after 
    % making sure that the mcmcm_simbio is in the MATLAB path. 
    %  Setup the paths to the various subdirectories of the toolbox, in particular the 
    % projects forlder, the source code, txtl modeling toolbox, and GWMCMC, etc. 

% get the path of this file
fp = mfilename('fullpath');

slashes = regexp(fp, '/');
% looks like slashes =
%      1     7    20    28    38    49    54    64    76

mcmc_simbio_relpath = fp(1:slashes(end)-1); 
% looks something like /Users/vipulsinghal/Dropbox/Documents/vipul_repo/mcmc/code_mcmc/mcmc_simbio

addpath([mcmc_simbio_relpath '/models_and_supporting_files']);
addpath([mcmc_simbio_relpath '/projects']);
addpath([mcmc_simbio_relpath '/exp_data']);


% need genpath only for the gwmcmc. the txtl modeling toolbox stored there
% is the default toolbox, but is not necessarily the one used. 
addpath([mcmc_simbio_relpath '/src']); 
addpath([mcmc_simbio_relpath '/src/gwmcmc']); 
end

