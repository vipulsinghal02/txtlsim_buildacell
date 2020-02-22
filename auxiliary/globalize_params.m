function [m] = globalize_params(m)
% globalize_params change the scope of all the parameters from the reaction
% to the model.
%   delete the locally scoped parameters and add parameters at the model
%   object level.
% Vipul Singhal, Nov 2017

% main loop over all the reactions
R = length(m.reactions);
for r = 1:R
    
    rx = m.reactions(r); % reaction class object
    
    % check if the reaction has an internal parameter object.
    if ~isempty(rx.KineticLaw.Parameters) 
        % nothing to do if there are no reaction scoped parameters.
        
        % get the parameter names
        %     parnames = get(rx.kineticlaw.getparameters, 'Name'); % dont need to use this.
        %     use this instead:
        parnames = rx.KineticLaw.ParameterVariableNames;
        
        % for each paramter name, get its value, then delete it.
        P = length(parnames);
        pvals = zeros(1, P);
        for p = 1:P
            % get the parameter values (to copy over to the global params)
            pvals(p) = rx.KineticLaw.Parameters(p).Value;
            % alt: get(rx.kineticlaw.getparameters, 'Value'), this returns a
            % cell array of numbers, also put this outside the P loop.
        end
        
        for p = 1:P
            % delete the parameter objects
            pTarget = sbioselect(rx, 'type', 'parameter', 'Name', parnames{p});
            pTarget.delete;
        end
        % Note that deleting the parameter objects does not remove the
        % parameter variable name property of the KineticLaw object.
        
        % now if reversible, setup the Kd parameter and the rule that Kd=kr/kf
        if rx.Reversible
            
            parname_base = parnames{1}(1:end-2);

            % assume a form like TXTL_UTR_UTR1_F, and so we want to remove the '_F' at the end.
            %todo: check if this is the format of all the parameters. 
            if isempty(sbioselect(m,'Type','Parameter', 'Name', [parname_base '_Kd']))
                addparameter(m, [parname_base '_Kd'], pvals(2)/pvals(1)) ;
            end
            
            if isempty(sbioselect(m,'Type','Parameter', 'Name', parnames{1}))
                addparameter(m, parnames{1}, pvals(1));
            end
            
            if isempty(sbioselect(m,'Type','Parameter', 'Name', parnames{2}))
                addparameter(m, parnames{2}, pvals(2)) ;
            end
            
            ruleStr = [parnames{2} ' = ' parname_base '_Kd*' parnames{1}];
            if isempty(sbioselect(m,'Type','Rule', 'Rule', ruleStr))
                addrule(m, ruleStr, 'initialAssignment');
            end
        else
            % irreversible reaction, just add the parameter to the model scope.
            %
            if isempty(sbioselect(m,'Type','Parameter', 'Name', parnames{1}))
                addparameter(m, parnames{1}, pvals(1));
            end
            
        end
    end
    
end

end

