classdef txtl_reaction_config
    %TXTL_REACTION_CONFIG Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        NTPmodel;
        AAmodel;
        Transcription_Rate;
        Translation_Rate;
        DNA_RecBCD_Forward;
        DNA_RecBCD_Reverse;
        DNA_RecBCD_complex_deg;
        Protein_ClpXP_Forward;
        Protein_ClpXP_Reverse;
        Protein_ClpXP_complex_deg;
        RNAP_S70_F;
        RNAP_S70_R;
        GamS_RecBCD_F;
        GamS_RecBCD_R;
        TL_AA_Forward;
        TL_AA_Reverse;
        TL_AGTP_Forward;
        TL_AGTP_Reverse;
        Ribosome_Binding_F;
        Ribosome_Binding_R;
        RNA_deg;
        RNase_F;
        RNase_R;
        NTP_Forward_1;
        NTP_Reverse_1;
        NTP_Forward_2;
        NTP_Reverse_2;
        AGTP_Concentration;
        CUTP_Concentration;
        AA_Concentration;
        RNAPbound_termination_rate;
        Ribobound_termination_rate;
        RecBCD_ic;
        RNAP_ic;
        Ribo_ic;
        RNase_ic;
        ATP_degradation_rate;
        ATP_regeneration_rate;
        ATP_degradation_start_time
    end
    
    methods
        function rConf = txtl_reaction_config(name)
            listOfProperties = properties(rConf);
            fileName = [name '_config.csv'];
            if exist(fileName,'file') == 2
                fid = fopen(fileName, 'rt');
                % current parameter file format
                % Param_name, Param_type, Value, Comment
                M =  textscan(fid,'%q%q%q%q','EndOfLine','\r\n', 'delimiter',',','CollectOutput',true);
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

