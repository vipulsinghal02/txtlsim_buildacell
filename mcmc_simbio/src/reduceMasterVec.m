function [rmv, rpr, sgnames] = reduceMasterVec(master_info)
	% in mcmc_info = mcmc_info_constgfp3ii(modelObj), we have the note: 
	% semanticGroups = {1, [2 4] [3 5]}; % cant do this, then the points never
% get differentiated at all. need some jitter. think about this actually.
% 
% return to that and try some things out. % working on this 
% nov 2018
% Copyright, Vipul Singhal
	mv = master_info.masterVector;
	estParamsIx = setdiff((1:length(mv))', master_info.fixedParams);
	logp = mv(estParamsIx);
	% reduce the logp (take only the first element for each group)
    % nreduc = sum(cellfun(@numel, master_info.semanticGroups)) ...
    %             - numel(master_info.semanticGroups);
    % reducedLength = length(logp) - nreduc;
    reducedLength = length(master_info.semanticGroups);
    % reduced master vector (actually the to-be-estimated part of it)
    rmv = zeros(reducedLength, 1);
    % reduced parameter ranges
    pr = master_info.paramRanges;
    rpr = zeros(reducedLength, 2);
    sgnames = cell(reducedLength, 1); % list of names for the sematic groups
    mastnames = master_info.estNames;
	for i = 1:reducedLength 
		sgi = master_info.semanticGroups{i}; % ith sematic group
		rmv(i) = logp(sgi(1));
		rpr(i, :) = [max(pr(sgi, 1)) min(pr(sgi, 2))];
		sgnames{i} = mastnames{sgi(1)}; % the first name in the group
	end

end



