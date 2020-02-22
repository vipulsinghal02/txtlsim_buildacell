function [da, mi, tv] = exportmobj(mi, data_info, fitOption)
%exportmobj export simbiology model object and update the model info struct
%with the ordering indices and the exported model object. 

nTopo = length(mi);
nGeom = zeros(length(mi),1);
for i = 1:nTopo
    nGeom(i) = length(mi(i).dataToMapTo);
end

% V2 data
% for each topology geometry pair, compute the data to fit
% ASSUME that the data dimensions across geometries is the same
% only differs across topologies at most. 
% Despite this, we allow for different geometries within a topology 
% to point to different data info array elements, just that all of 
% these elements must have the same number of timepoints, measured species, 
% replicates and dosing combinations. (version 3 of this code can be even 
% more general, with each topo-geom pair getting its own cell. this will be 
% slower.)

da = cell(nTopo, 1);
emo = cell(nTopo, 1);
tv = cell(nTopo, 1);
for i = 1:nTopo % each topology
    % each of the nGeom(i) geometries has a data info array element it points to. 
    % since the dimensions of these is assumed to be equal, we just use the first
    % one to set the empty array: 
    data_info_element = mi(i).dataToMapTo(1);
    currda = data_info(data_info_element).dataArray;
    % Transform Experimental - compute mean or median or nothing
    currda = computeFitOption(currda,  fitOption);
    da{i} = currda;
    tv{i} = data_info(data_info_element).timeVector;
    for j = 2:nGeom(i) % each geometry
        data_info_element = mi(i).dataToMapTo(j);
        currda = data_info(data_info_element).dataArray;
        % Transform Experimental - compute mean or median or nothing
        currda = computeFitOption(currda,  fitOption);
        % concatenate in the 5th dimension (the geometries dimension.)
        da{i} = cat(5, da{i}, currda); 
    end
    % EXPORT MODEL object to get it ready for MCMC
    % the resulting object is of class SimBiology.export.Model
    % documentation: 
    % https://www.mathworks.com/help/simbio/ref/simbiology.export.model-class.html
    % Sven Mesecke's blog post on using the exported model class for 
    % parameter inference applicaton. 
    % http://sveme.org/how-to-use-global-optimization-toolbox-algorithms-for-
    % simbiology-parameter-estimation-in-parallel-part-i.html

    
    mobj = mi(i).modelObj;

    enuo = mi(i).namesUnord; % estimated names unordered
    
    ep = sbioselect(mobj, 'Type', 'parameter', 'Name', ...
        mi(i).namesUnord);% est parameters

    es = sbioselect(mobj, 'Type', 'species', 'Name', ...
        mi(i).namesUnord);% est species
    
    
    aps = [ep; es]; % active parameters and species

    % reorder the parameter and species so they are in the same order as that
    % in the model. 
    eno = cell(length(aps), 1);% est names ordered
    ds = sbioselect(mobj, 'Type', 'species', 'Name', mi(i).dosedNames);
    
    % bugfix -- APRIL 19 2019: need to also reorder the dose values to
    % match the ordering in the model objects value info. See for example: 
    %     
    % ds
    % 
    %    SimBiology Species Array
    % 
    %    Index:    Compartment:    Name:                    InitialAmount:    InitialAmountUnits:
    %    1         contents        DNA ptet--utr1--deGFP    0                 
    %    2         contents        aTc                      0                 
    %    3         contents        DNA plac--utr1--tetR     0                 
    % 
    % mi(3).dosedNames
    % 
    % ans =
    % 
    %   3×1 cell array
    % 
    %     {'DNA ptet--utr1--deGFP'}
    %     {'DNA plac--utr1--tetR' }
    %     {'aTc'                  }
    
    ndoses = length(mi(i).dosedNames);
    doseOrderingIx = zeros(ndoses, 1);
    for q = 1:ndoses
        
        doseOrderingIx(q) = find(cellfun(@(x) strcmp({ds(q).Name}, x),...
            mi(i).dosedNames,...
            'UniformOutput', true));
        %IE, mi(i).dosedNames(doseOrderingIx(q)) == ds(q).Name
        % ie, use doseOrderingIx to transform mi(i).dosedNames ---> the
        % ordering in ds(q).Name, ie, expected by the exported model
        % object. 
    end
    
    
    emo{i} = export(mobj, [ep; es; ds]); % exported model object, dosed species names. 
    SI = emo{i}.SimulationOptions;

    % each of the nGeom(i) geometries has a data info array element it points to. 
    % since the dimensions of these is assumed to be equal, we just use the first
    % one to set the empty array: 
    data_info_element = mi(i).dataToMapTo(1);
    SI.StopTime = data_info(data_info_element).timeVector(end);
    accelerate(emo{i});
    
    mi(i).emo = emo{i}; % exported model object. 
    orderingIx = zeros(length(aps),1);
    orderingIx2 = orderingIx;
    for k = 1:length(aps)
        eno{k} = aps(k).Name;
        for kk = 1:length(enuo)
            if strcmp(eno{k}, enuo{kk} )
                orderingIx(k) = kk; % eno = enuo(orderingIx);
                % the kth element of orderingIx is kk. so the kth element of 
                % enuo(orderingIx) is enuo(kk). But this is just eno(k). And eno 
                % has the property of the kth element being eno(k). (as seen 
                % from "if eno{k} == enuo{kk} ")

                orderingIx2(kk) = k; %i.e., enuo = eno(orderingIx2); 
                % the kkth element of orderingIx2 is k. so the kk th element of
                % eno(orderingIx2) is eno(k). But the vector with this property is 
                % simply enuo. (as seen from "if eno{k} == enuo{kk} ")
            end
        end
    end

    mi(i).orderingIx = orderingIx; % these two arrays will be VERY useful. 
    mi(i).orderingIx2 = orderingIx2; % this one being the second. 
    mi(i).namesOrd = eno; % est names ordered. 
    mi(i).doseReordering = doseOrderingIx; % dose names that were used to reorder the doses. 
    mi(i).dosedNames = mi(i).dosedNames(doseOrderingIx);
    mi(i).dosedVals = mi(i).dosedVals(doseOrderingIx,:);
end
end

