function listmodelparts(m)
%listmodelparts list the insides of a model
%   list the species, the reactions, the global parameters, the parameters
%   and which reactions they are tied to, the rules and the events. 
% input: model object, m. 

m.species
m.Reactions
m.Parameters
gp = getparam(m);
gp(:, [1, 3])
m.Rules
m.Events
end

