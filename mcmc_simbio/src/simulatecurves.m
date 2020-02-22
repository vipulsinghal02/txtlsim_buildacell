function [da, idxnotused] = simulatecurves(em,m, nSimCurves, dose, tv, ms, varargin)

p = inputParser;
p.addParameter('sbioplot', false, @islogical);
p.parse(varargin{:})
p = p.Results;

nICs = size(dose, 1);
nMS = length(ms);
kknotused = zeros(nSimCurves, 1);
da = nan(length(tv), nMS, nSimCurves, nICs);
for i = 1:nICs
    for kk=1:nSimCurves
        try
            sd = simulate(em, [exp(m(end - kk+1,:)'); dose(i,:)']);
            sd = resample(sd, tv);
            % since measured species is a cell array of cell arrays of strings,
            % we need to loop over the outer cell, summing the values of the
            % species in the inner cell.
            for ss = 1:nMS
                % ms{ss} is a cell array of strings
                [~, XX] = selectbyname(sd, ms{ss});
                X = sum(XX, 2);
                % output the relevant data to the samples
                da(:,ss,kk, i) = X;
            end
        catch ME
            ME.message
            kknotused(kk) = 1;
        end
    end
    
end

idxnotused = find(kknotused);
if p.sbioplot
    sbioplot(sd)
end

end

