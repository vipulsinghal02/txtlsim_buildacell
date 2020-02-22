function [mv] = rebuildMasterVec(rmv, mai)
	% rmv = reduced master vector or array of vectors
	% mai: master_info struct
	% mv = rebuilt master vector. 

	% sg = cell array of semantic group vectors
	sg = mai.semanticGroups;
	mv = zeros(size(mai.paramRanges, 1), size(rmv, 2)); 
	% size(rmv, 2) should be 1 right? No! This can also be 
	% of lenght number of walkers! so for example, minit can be 
	% fed into this function, as done by the integrableLHS_v2.m

	for i = 1:length(sg)
        % build the matrix of semanticGroup normalizations
        % This is very interesting! can write a paper about this. 
        % data initialization structure. 
        % This is a cool bit of code! Very non intuitive why I am doing
        % this here. Worth explaining somewhere for future reference.... 
        % but am I going to?
        % 
        %
        %
%         y = datasample(s,1:2,2,'Replace',false)
        % 
        multfact = ones(numel(sg{i}), size(rmv, 2));
        for k = 1:size(rmv, 2)
            y = datasample(1:numel(sg{i}),numel(sg{i}),'Replace',false);
            for j = 1:length(y)
                % work on the y(j)th index
                multfact(y(j), k) = (-0.25 + 0.5*rand(1) + 0.74)^(y(j)-1);
                
            end 
        end
        
		mv(sg{i}, :) = multfact.*repmat(rmv(i, :), numel(sg{i}), 1);
	end
end
