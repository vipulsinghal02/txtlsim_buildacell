function [m, varargout] = set_mobj_from_model_info(m, mi, mai, geomID)
%set_mobj_from_model_info : Set the parameters in an unexported model object
%from the parameters specified in a particular geometry. (assume the
%model_info has been pre-sliced to the desired topology.)
% Note that the parameters in the master_info struct are expected to be log
% transformed.
%
%

pnames = mi.namesUnord;
mv = mai.masterVector;
pvals = mv(mi.paramMaps(:, geomID));

% set parameter values in the model
paramsNotSet = [];
for i = 1:length(pnames)
    obj = sbioselect(m, 'name', pnames{i});
    if ~isempty(obj)
        switch class(obj)
            case 'SimBiology.Parameter'
                obj.Value = exp(pvals(i));
            case 'SimBiology.Species'
                obj.InitialAmount = exp(pvals(i));
        end
    else
        paramsNotSet = [paramsNotSet; pnames(i)];
        
        warning(['Object with name ' pnames{i} ' not found. Skipping setting its value.'])
    end
end
if nargout ==2
    varargout{1} = paramsNotSet;
end

end

