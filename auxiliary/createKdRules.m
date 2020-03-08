function createKdRules(Mobj, rxstr, pname, Kdval, Fval)
% This function takes a mobj, a reversible reaction string, deletes all
% parameters associated with the reaction, and uses another ordered list
% of parameter names, Kd values and F rate values to make globally scoped
% parameters associated with those reactions. it adds a rule which ties the
% F and R rates via the Kd 
% 
% m: model object
% rxstr: a string for the exact reaction whose parameters are are trying to bind
% pname: a cell of string that is an identifier for the parameters 
% Kdval and Fval are the Kd and forward rate parameter, used to determine
% the reverse rate parameter


rx = sbioselect(Mobj, 'reaction', ...
    rxstr); %
parnames = get(rx.kineticlaw.getparameters, 'Name');
for i = 1:length(parnames)
pTarget = sbioselect(rx, 'type', 'parameter', 'Name', parnames{i});
pTarget.delete;
end
rx.KineticLaw.ParameterVariableNames = {[pname '_F'],[pname '_R']};

if isempty(sbioselect(Mobj,'Type','Parameter', 'Name', [pname '_Kd']))
addparameter(Mobj, [pname '_Kd'], Kdval) ;
end

if isempty(sbioselect(Mobj,'Type','Parameter', 'Name', [pname '_F']))
addparameter(Mobj, [pname '_F'], Fval);
end

if isempty(sbioselect(Mobj,'Type','Parameter', 'Name', [pname '_R']))
addparameter(Mobj, [pname '_R'], 0) ;
end
ruleStr = [pname '_R = ' pname '_Kd*' pname '_F'];
if isempty(sbioselect(Mobj,'Type','Rule', 'Rule', ruleStr))
	addrule(Mobj, ruleStr, 'initialAssignment');
end
end

