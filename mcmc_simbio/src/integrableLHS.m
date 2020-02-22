function int_minit = integrableLHS(eMO, nW, paramranges,...
    enames, ds, varargin)
%integrableLHS generate a set of integrable latin hypercube distributed
% parameter points for simbiology models
% eMO = exported model object
% nW = number of walkers

% OLD VERSION:
% spread = log-spread of the parameter values around logp
% logp = the parameters that set the center of the latin hypercube
% NEW VERSION:
% just specify the parameter ranges explicitly. these are log transformed. 

% enames = names of estimated parameters
% ds = dosing strategy
% optional name value pair: 'multiopt_params', mop is a matrix containing
% the indiced for the parameters to use for each sub optimization problem.
% it has rows corresponding to each sub optimization problem, and each row
% us the indices of the parameters to be estimated for that problem,
% followed by zeros to pad.


p = inputParser;
p.addParameter('multiopt_params', [], @isnumeric)
p.addParameter('distribution', 'LHS', @ischar)
p.addParameter('width', .001, @isnumeric)

p.parse(varargin{:});

p = p.Results;

if iscell(eMO)
    % TODO: update this multiopt version. for now, we just do the single
    % opt version. 
    mop = p.multiopt_params;
    assert(~isempty(mop),...
        'Please ensure that in multi optimization mode, the estimation structure is specified')
    assert(isequal(length(logp), length(estNamesFull), size(mop,2),...
        'The length of logp, full estimation names, and the # columns in estimation structure must match'));
    assert(isequal(length(eMO), size(mop, 1)), ...
        'The number of sub optimization problems must be consistent in the array of model objects and the estimation structure array');
    nOpt = length(eMO);
    nparam = length(logp); % logp not defined?
    % for each sub problem, get indices of the integrable points
    % Most of the time we seem for have 70% or more integrability. So assume we lose 30% of points per sub problem.
    npts = round(nW/max([0.65^nOpt 0.1])); %we allow for up to 10X multiplication, need to rethink the method if this does not suffice
    lhsamp = paramranges*(lhsdesign(npts, nparam)-0.5);
    minit=bsxfun(@plus,logp,lhsamp');
    IP = cell(nOpt,1);
    % number of points to generate initially should be
    tic
    for kk = 1:nOpt
        % here we reorder the parameters for the individual optimization problem. The estNames input to the sbiointegrable function below
        % should be in the order that the parameters appear in the exported model object.
        [~,~,V] = find(mop(kk,:));
        en = enames(V);
        % order en correctly by working through valueInfo
        count = 1;
        for ii = 1:length(eMO.ValueInfo)
            index1 = strcmp(en, eMO.ValueInfo(ii).Name);
            index2 = find(index1);
            if index2~= 0
                en2{count} = en{index2};
                count = count+1;
            end
        end
        
        IP{kk} = sbiointegrable(eMO{kk}, minit, en2, ds{kk});
        numInt{kk} = sum(all(IP{kk}==1, 2));
        intIx{kk} = find(all(IP{kk}==1, 2));
    end
    toc ;

    simul_int = (1:npts)';
    for kk = 1:nOpt
        simul_int = intersection(intIx{kk}, simul_int);
    end
    
    if nW >= length(simulint)
        int_minit = minit(:, simul_int);
        warning('Number of integrable points less than number specified. Using only integrable points.')
        
    else
        % randomly generate nW samples from list of integrable indices
        r = randperm(length(simulint), nW);
        int_minit = minit(:, simul_int(r));
    end
    
else
    
    nparam = length(enames);
    npts = round(nW*2); % can tolerate up to 50% non integrability.
    % Increase the factor here if you need to tolerate more.
    
    % generate latin hyper cube distributed points to test for
    % integrability
    switch p.distribution
        case 'LHS'
            lhsamp = lhsdesign(npts, nparam);
            lhsamp = lhsamp'; % nparam by npts matrix of LHS points
            
            minit= ...
                lhsamp.*(repmat(paramranges(:, 2), 1, npts)-repmat(paramranges(:, 1), 1, npts))+...
                repmat(paramranges(:, 1), 1, npts);
            
        case 'gaussian'
            midpt = (repmat(paramranges(:, 2), 1, npts) +...
                repmat(paramranges(:, 1), 1, npts))/2;
            
            minit = p.width*randn(nparam,npts)-p.width/2 + midpt;
                
        case 'unifrand'
            midpt = (repmat(paramranges(:, 2), 1, npts) +...
                repmat(paramranges(:, 1), 1, npts))/2;
            width = (repmat(paramranges(:, 2), 1, npts) -...
                repmat(paramranges(:, 1), 1, npts));
            minit = width.*rand(nparam,npts)-width/2 + midpt;
    end
    
    
   tic
    IP = sbiointegrable(eMO, minit, enames, ds);
    toc ;
    
    numInt = sum(all(IP==1, 2));
    intIx = find(all(IP==1, 2));
    if nW >= numInt
        int_minit = minit(:, intIx);
    else
        % randomly generate nW samples from list of integrable indices
        r = randperm(numInt, nW);
        int_minit = minit(:, intIx(r));
    end
end

end

