function createFrateRules(Mobj, rxstr, varargin)
% This function takes a mobj, a reaction string, deleted the parameter associated
% with the reaction, and creates a global parameter that informs the rate fo that reaction
% 
% m: model object
% rxstr: a string for the exact reaction whose parameters are are trying to bind
% pname: a string that has the letters F, R or Kd appended to it, and is an
% identifier for the parameters 
rx = sbioselect(Mobj, 'reaction', ...
    rxstr); %
	
if rx.reversible 
	error('Expecting an irreversible reaction')
else
	parname = get(rx.kineticlaw.getparameters, 'Name');
	Fval = get(rx.kineticlaw.getparameters, 'Value');
end

pTarget = sbioselect(rx, 'type', 'parameter', 'Name', parname);
pTarget.delete;

nvararg = length(varargin);
optargs = {parname, Fval};
optargs(1:nvararg) = varargin;
[pname, Fval] = optargs{:};

% when there is only one parameter in the given reaction, parnames is a string

rx.KineticLaw.ParameterVariableNames = pname; 
% I am guessing this needs to be a string. 

	if isempty(sbioselect(Mobj,'Type','Parameter', 'Name', pname))
		addparameter(Mobj, pname, Fval);
	end
end