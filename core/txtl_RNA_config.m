classdef txtl_RNA_config
    %TXTL_REACTION_CONFIG Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        att1_anti1_termination_rate;
        att1_anti1_repression_F_rate;
        att1_anti1_repression_R_rate;
        att1_auto_termination_rate;
        att1_anti2_Xtalk_F_rate;
        att1_anti2_Xtalk_R_rate;
        att1_anti2_Xtalk_termination_rate;  
        att2_anti2_termination_rate;
        att2_anti2_repression_F_rate;
        att2_anti2_repression_R_rate;
        att2_auto_termination_rate;
        att2_anti1_Xtalk_F_rate;
        att2_anti1_Xtalk_R_rate;
        att2_anti1_Xtalk_termination_rate; 
        att1_anti11_repression_F_rate;
        att1_anti11_repression_R_rate;
        att1_anti11_termination_rate;
        att2_anti1_anti1_Xtalk_F_rate;
        att2_anti1_anti1_Xtalk_R_rate;
        att2_anti11_Xtalk_termination_rate;
    end
    
    methods
        function rConf = txtl_RNA_config(name)
            listOfProperties = properties(rConf);
            fileName = [name '_config.csv'];
            if exist(fileName,'file') == 2
                fid = fopen(fileName);
                % current parameter file format
                % Param_name, Param_type, Value, Comment
                M =  textscan(fid,'%q%q%q%q', 'delimiter',',','CollectOutput',true);
                M = M{1};
                fclose(fid);
                
                % set the object properties from the file.
                % Numeric/Expression type parameters are distinguished.
                for k = 1:size(listOfProperties)
                    index = find(cellfun(@(x) strcmp(x,listOfProperties(k)),M(:,1)) > 0);
                    if index > 0
                        if strcmp(M{index,2},'Numeric')
                            eval(sprintf('rConf.%s = %s;',M{index,1},M{index,3}));
                        else
                            % trying to evaluate expressions, if fails
                            % expression saved as string for later
                            % clarification
                            try
                               eval(sprintf('rConf.%s = %s;',M{index,1},M{index,3})); 
                            catch err
                               eval(sprintf('rConf.%s = ''%s'';',M{index,1},M{index,3}));
                            end
                        end
                    end
                end
            else
                error('the file: %s does not exist!',name);
            end
        end
        
        
      
    end
    
end

