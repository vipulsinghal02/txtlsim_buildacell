function int_minit = integrableLHS_v2(mi, mai, ri, varargin)
    % version 2 of the integrable LHS function. 
    % 
    % OLD HELP FILE:
    %integrableLHS generate a set of integrable latin hypercube distributed
    % parameter points for simbiology models
    % eMO = exported model object
    % ri.nW = number of walkers

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

    p.addParameter('Parallel',true,@islogical);
    p.parse(varargin{:});

    p = p.Results;

    % compute the reduced number of parameters.
    % nreduc = sum(cellfun(@numel, mi.semanticGroups)) ...
    %             - numel(mi.semanticGroups);

    % reducedvec -> mastervec -> -> -> distribute across all models 
    % and geometries and doses. if the integration passes for all cases, 
    % then that points is integrable. Probably need to expand the candidate 
    % points to be 5x the original number of walkers. 
    % make sure the parameter order is correct. 
    % 
    [~, rpr, ~] = reduceMasterVec(mai);
    nparam = size(rpr, 1);


    npts = round(ri.nW*3); % can tolerate up to 67% non integrability.
    % Increase the factor here if you need to tolerate more.
    
    % Compute the parameter sharing across all topologies and geometries for 
    % initial walker estimation purposess. 

    % generate latin hyper cube distributed points to test for
    % integrability

    switch p.distribution
        case 'LHS'
            lhsamp = lhsdesign(npts, nparam);
            lhsamp = lhsamp'; % nparam by npts matrix of LHS points
            % Note that this is a npts dimensional hypercube, with sides
            % from 0 to 1, excluding the edges. 
            
            
            minit= ...
                lhsamp.*(repmat(rpr(:, 2), 1, npts)-...
                    repmat(rpr(:, 1), 1, npts))+...
                repmat(rpr(:, 1), 1, npts);
            
        case 'gaussian'
            midpt = (repmat(rpr(:, 2), 1, npts) +...
                repmat(rpr(:, 1), 1, npts))/2;
            
            minit = p.width*randn(nparam,npts)-p.width/2 + midpt;
                
        case 'unifrand'
            midpt = (repmat(rpr(:, 2), 1, npts) +...
                repmat(rpr(:, 1), 1, npts))/2;
            width = (repmat(rpr(:, 2), 1, npts) -...
                repmat(rpr(:, 1), 1, npts));
            minit = width.*rand(nparam,npts)-width/2 + midpt;
    end
    
    % rebuild the master vector
    minit_justEstParams = rebuildMasterVec(minit, mai); 
    % only has #estparams in the columns. 
    % not the full master vector. 
    % all the values should be in log space. 

    % we build the full set of npts master vectors, arranged into a matrix
    % of size #full master vector elements x npts: 
    mv = mai.masterVector;
    minit_fixedAndEstParams = repmat(mv, 1, npts);
    estParamsIx = setdiff((1:length(mv))', mai.fixedParams);
    minit_fixedAndEstParams(estParamsIx, :) = minit_justEstParams; 
    % i changed this line on 3.12.2018. not sure if this break previous
    % code or fixes it. check how the previous code was working in the
    % first place. yeah i think the previous code, which was only the vnprl
    % and the protein with the mrna parameters fixed to 10 sets of values
    % all never had semantic groups for parameters to be ESTIMATED. 

    % each column of minit, when distributed across topologies and geometries
    % should work for every topology geometry pair. 

    IProw_old = ones(1, npts);
    for i = 1:length(mi) % for each topology

        %number of params in a given topo (model)
        nParam_TopoGeom = size(mi(i).paramMaps, 1); 

        % define dose matrix. This is correct, can also put a {} around the dosedvals 
        % matrix. Cant remove the braces around the names. 
        ds = struct('names', {mi(i).dosedNames}, 'dosematrix', mi(i).dosedVals);


        for j = 1:size(mi(i).paramMaps, 2) % for each geometry, do everything. 
            
            % build the minit for this topo - geom pair
            % from the master one (ie, minit_fixedAndEstParams)
            pIX_TopoGeom = mi(i).paramMaps(mi(i).orderingIx, j);
            minitTopoGeom = minit_fixedAndEstParams(pIX_TopoGeom, :);

            % REORDER TO MAKE THEM CORRECT WITH THE EXPORTED MODEL'S
            % parameter ORDERING. 
            % 
            % ??? MAYBE NOT NEEDED. JUST REORDER AFTR THE ESTIMATION. THE
            % PARAM RANGERS ARE MAYBE ALREADY ORDERED. 
            
%             minitTopoGeom_reordered = minitTopoGeom(
            
            
            
            % now simulate this the t-g for this minit and
            % for each column of minit, report if it passes.
            % IP is a matrix of dimension npts x # dose combinations. 
            % Three possible values: 
            % 0 = integration tol not met 
            % 1 = run successfully
            % 2 = unknown error. 
            % 
            tic
            IP = sbiointegrable(mi(i).emo, minitTopoGeom, mi(i).namesOrd, ds, p.Parallel);
            toc
            IP = IP'; % make IP nDoses x npts 


            % numInt = sum(all(IP==1, 1));
            IProw_new = all(IP==1, 1);
            % intIx = find(all(IP==1, 1));
            IPtemp = [IProw_old; IProw_new];
            IProw_old = all(IPtemp==1, 1);
        end
        % do this for each geometry and each topology. the set of points that
        % pass every t-g pair are our valid starting point. could be stringent. 
    end
    IProw = IProw_old;
    numInt = sum(IProw==1);
    disp([num2str(numInt) ' points out of ' num2str(size(minitTopoGeom,2)) ...
        ' are integrable. Need ' num2str(ri.nW) ' walkers.'])
    intIx = find(IProw==1);
    if ri.nW >= numInt
        % not enough integrable points
        warning(['Not enough integrable points. '...
            'Reducing to maximal set of integrable points.'])
        int_minit = minit_justEstParams(:, intIx);
    else
        % randomly generate ri.nW samples from list of integrable indices
        r = randperm(numInt, ri.nW);
        int_minit = minit_justEstParams(:, intIx(r));
    end


end

