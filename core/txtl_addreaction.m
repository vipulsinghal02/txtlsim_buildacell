function txtl_addreaction(tube,reactionEq,kineticLaw,parameters,varargin)

%%% Vesicule mode %%%
if tube.Userdata.Vesicule && size(tube.Compartments,1) > 1
    skipOperation = false;
    if nargin > 4 && ~strcmp(varargin{1},'membrane_transport')
        compID = varargin{1};
    elseif nargin <= 4
        compID = tube.Compartments(1).Name;
    else
        % the reaction equation is a membrane transport equation, all the
        % compartment names are provided -> no action needed
        skipOperation = true;
    end
    if ~skipOperation
        % spliting up the reatctionEQ
        [matchstr splitstr] = regexp(reactionEq,'+|<->|<-|->','match','split');
        % removeing spaces and adding compartment name
        splitstr = cellfun(@(x) [compID '.' strtrim(x)],splitstr,'UniformOutput',false);
        % adding space before and after operations (SimBio requirement)
        matchstr = cellfun(@(x) [' ' x ' '],matchstr,'UniformOutput',false );
        % spaceholder for reactionEq with Compartment id
        mergeStr = cell(1,size(splitstr,2) + size(matchstr,2));
        mergeStr(1:2:end) = splitstr;
        mergeStr(2:2:end) = matchstr;
        
        reactionEq = strjoin(mergeStr,'');
    end
    
end
%%% Vesicule mode %%%

Robj1 = addreaction(tube, reactionEq);
Kobj1 = addkineticlaw(Robj1, kineticLaw);

for k=1:size(parameters,1)
    addparameter(Kobj1, parameters{k,1}, parameters{k,2});
end

set(Kobj1, 'ParameterVariableNames', parameters(:,1)');


end