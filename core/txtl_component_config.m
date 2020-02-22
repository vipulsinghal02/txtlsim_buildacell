classdef txtl_component_config
    %TXTL_COMPONTENT_CONFIG Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        RNAPbound_Forward
        RNAPbound_Reverse
        RNAPbound_Forward_actv
        RNAPbound_Reverse_actv
        RNAPbound_Forward_2
        RNAPbound_Reverse_2
        RBS_Length
        spacer_Length
        Promoter_Length
        Thio_Length
        Junk_Length
        Gene_Length
        LVA_tag_Length
        Terminator_Length
        Dimmerization_Forward
        Dimmerization_Reverse
        Tetramerization_Forward
        Tetramerization_Reverse
        Protein_Inducer_Forward
        Protein_Inducer_Reverse
        Inducer_Degradation
        Protein_Inducer_Forward_2
        Protein_Inducer_Reverse_2
        Inducer_Degradation_2
        % N*2 matrix for possible promoter-protein complex interaction
        % we assume that there is 
        %   - N possible protein complex can bind to the promoter sites 
        %   - forward & reverse reaction rates
        DNA_Sequestration
        non_s70_factor_Forward
        non_s70_factor_Reverse
        Protein_Maturation
        activation_F
        activation_R
        generic_rate1 
        generic_rate2 
        generic_rate3
        %used for special reactions one may want to add that are only present 
        %in a couple of parts, like phosphorylation.
        prot_deg_F
        prot_deg_R
        prot_deg_unfold
        ClpX_deg
        Ribosome_Binding_F
        Ribosome_Binding_R
        Combinatorial_activ_knockoff_F
        Combinatorial_activ_knockoff_R
        
    end
    
    methods
        function compConf = txtl_component_config(name)
            listOfProperties = properties(compConf);
            fileName = ['txtl_param_' name '.csv'];
            if exist(fileName,'file') == 2
                fid = fopen(fileName, 'rt');
                % current parameter file format
                % Param_name, Param_type, Value, Comment
                M =  textscan(fid,'%q%q%q%q','EndOfLine','\r\n', 'delimiter',',','CollectOutput',true);
                M = M{1};
                fclose(fid);

                % set the object properties from the file.
                % Numeric/Expression type parameters are distinguished.
                for k = 1:size(M)
                    index = find(cell2mat(cellfun(@(x) strcmp(x,M(k,1)),listOfProperties,'UniformOutput',false)) > 0);
                    if index > 0
                        if strcmp(M{k,2},'Numeric')
                            eval(sprintf('compConf.%s = %s;',M{k,1},M{k,3}));
                        else
                            % trying to evaluate expressions, if fails
                            % expression saved as string for later
                            % clarification
                            try
                               eval(sprintf('compConf.%s = %s;',M{k,1},M{k,3})); 
                            catch err
                               eval(sprintf('compConf.%s = ''%s'';',M{k,1},M{k,3}));
                            end
                        end
                    else
                        r = regexp(M{k,1},'^DNA_Sequestration_C([1-9])+_(F|R)$','tokens','once');
                        oneComplex = regexp(M{k,1},'^DNA_Sequestration_(F|R)$','tokens','once');
                        if size(r,2) > 0
                            p = 1;
                            if strcmp(r{2},'R')
                                p = 2;
                            end
                            compConf.DNA_Sequestration(str2double(r{1}),p) = str2double(M{k,3});
                        elseif ~isempty(oneComplex)
                            if ~isempty(oneComplex{1}) && ischar(oneComplex{1})
                            p = 1;
                                if strcmp(oneComplex{1},'R')
                                    p = 2;
                                end
                            end
                            compConf.DNA_Sequestration(1,p) = str2double(M{k,3});
                        end 
                    end % end of if index > 0 
                end % end of for
            else
                error('the file: %s does not exist!',name);
            end
        end % end of constructor
        
        % possible inputs
        % - nth protein complex forward rates
        % getDNASequestrationRates(n,'Forward')
        % getDNASequestrationRates(n,'F')
        % - first protein complex (assumed), forward rates
        % getDNASequestrationRates('Forward')
        % getDNASequestrationRates('F')
        %
        function rate = getDNASequestrationRates(varargin)
            compConf = varargin{1};
            switch nargin
                case 2
                    Conf = 1;
                    direction = varargin{2};
                case 3 
                    Conf = varargin{2};
                    direction = varargin{3};
                otherwise 
                    error('Number of argument should be either 2 or 3, not %s',nargin);
            end
            
            if ischar(direction)
                f = strcmp({'Forward','F'},direction);
                r = strcmp({'Reverse','R'},direction);
                if sum(f) > 0
                    direction = 1;
                elseif sum(r) > 0
                    direction = 2;
                end
            end
            rate = compConf.DNA_Sequestration(Conf,direction);
        end
        
    end
    
end

