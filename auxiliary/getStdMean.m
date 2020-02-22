function [data_mean,data_std] = getStdMean(varargin)

% single argument: either one matrix or cell arrray of matrices
if nargin == 1
    if iscell(varargin{1})
        opMode = 'MultiMatrix';
        multiData = varargin{1};
        
    else
        opMode = 'singleMatrix';
        data = varargin{1};
    end
    % multiple matrices are given
else
    opMode = 'MultiMatrix';
    multiData = varargin;
end

% select the right operation mode and calculate first and second modes
switch opMode
    case 'singleMatrix'
        
        data_mean = mean(data,2);
        if size(data,2) == 1
            data_std = zeros(size(data));
        else
            data_std  = std(data')';
        end
        
    case 'MultiMatrix'
        dataDim = cell2mat(cellfun(@size,multiData,'UniformOutput',false));
        
        assert(all(dataDim(:,2) == dataDim(1,2)),'the number of wells should be the same in each experiment');
        
        NumDataPoints = min(dataDim(:,1));
        % dataDim is a 3D array
        if size(dataDim,2) == 3
            numOfChannels = min(dataDim(:,3));
        else
        % dataDim is a 2D array   
            numOfChannels = 1;
        end

        for l = 1:numOfChannels
        
            for k= 1:dataDim(1,2)
                data           = cell2mat(cellfun(@(x) x(1:NumDataPoints,k,l),multiData','UniformOutput',false)) ;
                data_mean(:,k,l) = mean(data,2);
                data_std(:,k,l)  =  std(data')';
            end
        end
    otherwise
        error('unknown operation mode');
end

end