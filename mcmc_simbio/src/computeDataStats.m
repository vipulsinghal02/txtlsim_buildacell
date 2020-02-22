function [summst, spreadst] = computeDataStats(dataArray, dispmode)
% da: data array
% datasummary: mean, median or none
% dataspread: curves, std or none.
% rD: replicates dimension

switch dispmode
    case 'mean'
        summst = mean(dataArray, 3);
        spreadst = [];
    case 'median'
        % sum over the time dimension
        % tD MUST be 1 for this to work.. I tried being fully general,
        % but what is the point? It's unnecessarily difficult.
        % compute the indexes of the median (in terms of sum / integral)
        % curves over the replicates
        [ix, mdvals] = medianIndex(sum(dataArray, 1), 3);
        % again, rD MUST be 3.
        summst = medianReplicate(dataArray, ix);
        % spreadstatistic is the empty vector
        spreadst = [];
    case 'meanstd'
        summst = mean(dataArray, 3);
        spreadst = std(dataArray, 0, 3);
    case 'meancurves'
        summst = mean(dataArray, 3);
        spreadst = dataArray;
    case 'medianstd'
        [ix, mdvals] = medianIndex(sum(dataArray, 1), 3);
        summst = medianReplicate(dataArray, ix);
        spreadst = std(dataArray, 0, 3);
    case 'mediancurves'
        [ix, mdvals] = medianIndex(sum(dataArray, 1), 3);
        summst = medianReplicate(dataArray, ix);
        spreadst = allButMedianCurve(dataArray, ix);
    case 'curves'
        summst = [];
        spreadst = dataArray;
        
    otherwise
        error(['Invalid data display mode. Must be one of: ''mean'','...
            ' ''median'' , ''meanstd'',''medianstd'', ''meancurves'','...
            ' ''mediancurves'', ''curves''.'])
end
end
