function currda = computeFitOption(da, fo)
    % fo is the fit option. 
    % da is the data array
    switch fo
        case 'FitMedian'
            % Compute the curvewise median of the data. 
            [ix, mdvals] = medianIndex(sum(da, 1), 3); 
            currda = medianReplicate(da, ix);
        case 'FitMean'
            currda = mean(da, 3);
        case 'FitAll'
            currda = da;
        otherwise
            error(...
                ['Invalid fit option. Read the documentation'...
                ' for how to specify inputs.'])
    end

end