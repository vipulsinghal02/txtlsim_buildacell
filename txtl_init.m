function txtldir = txtl_init
fp = mfilename('fullpath'); 
slashes = regexp(fp, '/');
filedir = fp(1:slashes(end)-1);
addpath(filedir);
addpath(genpath([filedir '/auxiliary']))
addpath([filedir '/components'])
addpath([filedir '/config'])
addpath([filedir '/core'])
addpath([filedir '/examples'])
addpath([filedir '/examples/CompanionFiles'])
addpath([filedir '/tests'])
addpath([filedir '/data'])
addpath([filedir '/mcmc_simbio'])
txtldir = filedir;
end
