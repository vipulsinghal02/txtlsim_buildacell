function txtl_enzyme_resource_degradation(modelObj)
%


if isfield(modelObj.UserData, 'energymode') && strcmp(modelObj.UserData.energymode, 'regeneration')
        atp_deg_rate = modelObj.UserData.ReactionConfig.ATP_degradation_rate;
        atp_regen_time = modelObj.UserData.ReactionConfig.ATP_degradation_start_time;
        atp_reg_rate = modelObj.UserData.ReactionConfig.ATP_regeneration_rate;
        % After some time, ATP regeneration stops, leading to an overall decrease in
        % ATP concentrations. c.f. V Noireaux 2003
        parameterObj = addparameter(modelObj, 'AGTPreg_varying', atp_reg_rate, 'ConstantValue', false);
        
        % time of the regenration system turn off
        parameterObj = addparameter(modelObj, 'AGTPdeg_time', atp_regen_time, 'ConstantValue', true);
        
        parameterObj = addparameter(modelObj, 'AGTPreg_ON', atp_reg_rate, 'ConstantValue', true);
        parameterObj = addparameter(modelObj, 'AGTPdeg_rate', atp_deg_rate, 'ConstantValue', true);
        
        
        evt2 = addevent(modelObj, 'time <= AGTPdeg_time' , 'AGTPreg_varying = AGTPreg_ON');
        
        evt3 = addevent(modelObj, 'time > AGTPdeg_time',...
            'AGTPreg_varying = 0');
        
        % constantly active first order degradation. 
        reactionObj = addreaction(modelObj,'AGTP -> AGMP');
        kineticlawObj = addkineticlaw(reactionObj, 'MassAction');
        set(kineticlawObj, 'ParameterVariableName', 'AGTPdeg_rate');
        
        % regeneration system
        reactionObj = addreaction(modelObj,'AGMP -> AGTP');
        kineticlawObj = addkineticlaw(reactionObj, 'MassAction');
        set(kineticlawObj, 'ParameterVariableName', 'AGTPreg_varying');
        
else
    
    atp_deg_rate = modelObj.UserData.ReactionConfig.ATP_degradation_rate;
    atp_deg_time = modelObj.UserData.ReactionConfig.ATP_degradation_start_time;
    
    % After some time, ATP regeneration stops, leading to an overall decrease in
    % ATP concentrations. c.f. V Noireaux 2003.
    parameterObj = addparameter(modelObj, 'AGTPdeg_F', 0, 'ConstantValue', false);
    
    parameterObj = addparameter(modelObj, 'AGTPdeg_time', atp_deg_time, 'ConstantValue', true);
    
    parameterObj = addparameter(modelObj, 'AGTPdeg_rate', atp_deg_rate, 'ConstantValue', true);
    
    evt2 = addevent(modelObj, 'time <= AGTPdeg_time' , 'AGTPdeg_F = 0');
    
    evt3 = addevent(modelObj, 'time > AGTPdeg_time',...
        ['AGTPdeg_F = AGTPdeg_rate']);% '=' num2str(atp_deg_rate)]
    
    reactionObj = addreaction(modelObj,'AGTP -> AGTP_USED');
    
    kineticlawObj = addkineticlaw(reactionObj, 'MassAction');
    set(kineticlawObj, 'ParameterVariableName', 'AGTPdeg_F');
end

end