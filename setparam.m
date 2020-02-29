function p = setparam(m, rStr, pStr, pVal)
%setparam: Set parameters in a model object
%   m is model object, rStr is a single string of the reaction
% pStr is a single parameter name string, or a 2 by 1 cell of parameter
% strings
% pVal is a scalar, or a 2 by 1 vector of parameter values.

if iscell(pStr)
    if strcmp(rStr, 'global')
        param1 = pStr{1};
        param2 = pStr{2};
        p1 = sbioselect(m.parameters, 'Name',  param1);
        p2 = sbioselect(m.parameters, 'Name',  param2);
        set(p1, 'Value', pVal(1));
        set(p2, 'Value', pVal(2));
        p = [p1;p2];
    else
        param1 = pStr{1};
        param2 = pStr{2};
        rObj = sbioselect(m,'Type','reaction','reaction',rStr);
        p1 = sbioselect(rObj,'Type','parameter','Name',param1);
        p2 = sbioselect(rObj,'Type','parameter','Name',param2);
        set(p1, 'Value', pVal(1));
        set(p2, 'Value', pVal(2));
        p = [p1;p2];
    end
elseif ischar(pStr)
    % rStr needs to be a string, pVal needs to be a scalar
    if strcmp(rStr, 'global')
        param1 = pStr;
        p = sbioselect(m.parameters, 'Name',  param1);
        set(p, 'Value', pVal(1));
    else
        param1 = pStr;
        rObj = sbioselect(m,'Type','reaction','reaction',rStr);
        p = sbioselect(rObj,'Type','parameter','Name',param1);
        set(p, 'Value', pVal(1));
    end
else
    error('parameter strings must be cells or a string')
end


end

