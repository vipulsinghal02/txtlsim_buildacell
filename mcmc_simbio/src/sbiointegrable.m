function IP = sbiointegrable(m, pv, pn, ds, parllel)
%SBIOINTEGRABLE Check integrability of model on a list of specified parameter
%points
%
% m: is either a simbiology model object or an exported
% and accelerated version of a simbiology model object
%
% pv is an M x N matrix of log transformed parameter values, where M is the number of model
% parameters and N is the number of points in parameter space.
%
% pn is the list of parameters and species (initial conditions)
% corresponding to each column of pv.
%
% ds is the list of dosed species names and concentrations to use. This is
% a MATLAB struct. It is slightly different than the structure specified according to the fucntion
% txtl_gen_noiseless_data, however if a structure in that format is given as an input,
% this function is able to still parse that.
% The general format is ds = struct('names', {{'A';'B';'C'}}, 'dosematrix', dm)
%
%

% check inputs:
% 1) model can be exported or not
% 2) dosing strategy ds can be in the traditional format or in the matrix
% format.
if isa(m, 'SimBiology.Model')
    exported = false;
elseif isa(m, 'SimBiology.export.Model')
    exported = true;
else
    error('m must be a simbiology model')
end


% parllel = false; % in the future enable this as an argument 

if all(isfield(ds, {'species', 'concentrations'}))
%     dstype = 1; % type 1 dosing strategy is the usual type we use in the top level file
    % convert to type 2
    dnames = {};
    for i = 1:length(ds)
        dnames = [dnames; {ds(i).species}];
    end
    nds = length(ds); % number of species to dose
    ndc = length(ds(1).concentrations);% number of dose combinations
    
    dm = zeros(nds, ndc);
    for j = 1:ndc
        for i = 1:nds
            dm(i,j) = ds(i).concentrations(j);
        end
    end
    
ndc = length(ds(1).concentrations); % number of dose combinations

elseif all(isfield(ds, {'names', 'dosematrix'}))
%     dstype = 2;
    dnames = ds.names;
    dm = ds.dosematrix;
    [nds ndc] = size(dm);
    % dosematrix is a # of dose species x # of dose combinations matrix
    % its columns are used to augment the parameter value vector we
    % simulate the exported model with.
else
    error('Unknown dosing strategy struct format')
end


% error checking or exporting of the model
if exported
    % check if the pn and ds arrays match the exported model's internal
    % specification of what the parameters to vary and doses to use are.
    % first the parameters to be estimated, then the species initial
    % concentrations to be estimated, then the species initial
    % copncentrations to be dosed.
    fullpn = [pn; dnames];
    errmsg = ['the number of variables in the exported model object should'...
    ' be the sum of the number of parameters being varied and number of species to be dosed.'];
    assert(length(fullpn)==length(m.ValueInfo), errmsg)
    for i = 1:length(fullpn)
        errmsg2 = sprintf('user specified value %d is %s while model expects %s',...
            i, fullpn{i}, m.ValueInfo(i).Name);
        assert(strcmp(fullpn{i}, m.ValueInfo(i).Name), errmsg2)
    end
else
    % export the model, using pn and ds to specifiy the parameters and
    % species to keep variable.
    m1 = copyObj(m);
    clear m;
    m = export(m1, [pn; dnames]);
end

N = size(pv, 2); % number of parameter combinations (ie # of points in parameter space)


IP = nan(N, ndc); % Integrable Points

for d = 1:ndc % d is the index for the dose concentrations
    cdv = dm(:,d); % current vector of dose values
    % note: nov 26 2017: really dont need parallel computing here. we are just running this
    % once, right? Maybe not. I think there was a use where I was trying a
    % lot of cases: different step sizes, different atol, different rtols
    % etc. 
    if parllel == true
        disp(['Testing integrability for dose number ' num2str(d) '.']);
        parfor n = 1:N % n is the index for the number of total parameter combinations/points.
            cpv = pv(:,n); % current parameter values
            
            % recall that in the exported model object the format is
            % [parameters to simulate with ; species to s ; 
            fv = [exp(cpv); cdv] ;% full set of values to vary in the model for this simulation
            % run the model with the specified parameters
            try
                simulate(m, fv); % sd is simdata
                IP(n,d) = 1;
            catch ME
                if strcmp(ME.identifier, 'DESuite:ODE15S:IntegrationToleranceNotMet')
                    % store parameter vals
                    IP(n,d) = 0;
                else
                    % hopefully this never happens
                    warning('An unknown Error has occurred, Code 2')
                    IP(n,d) = 2; % Unknown error
                end
            end
        end
    else
        for n = 1:N % n is the index for the number of total parameter combinations.
            cpv = pv(:,n); % current parameter values
            fv = [exp(cpv); cdv] ;% full set of values to vary in the model for this simulation
            % run the model with the specified parameters
            try
                simulate(m, fv); % sd is simdata
                IP(n,d) = 1;
            catch ME
                if strcmp(ME.identifier, 'DESuite:ODE15S:IntegrationToleranceNotMet')
                    % store parameter vals
                    IP(n,d) = 0;
                else
                    % hopefully this never happens
                    warning('An unknown Error has occurred, Code 2')
                    IP(n,d) = 2; % Unknown error
                end
            end
        end
        
        
        
    end
    
    
    
    
    
end

