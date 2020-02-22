function [endtime_ind,endtime] = expression_endtime(t_vec,data,varargin)

% number of data channel in the data matrix
dataChannels = size(data,2);
%we have reached the total signal above the threshold level
threshold = 0.95;

for k=1:dataChannels
    if sum(data(:,k)) > 0
        indend = find(data(end,k)*threshold < data(:,k));
        endtime_ind(k) = indend(1);
        endtime(k) = t_vec(endtime_ind(k));
    else 
        endtime(k) = 0;
        endtime_ind(k) =0;
    end
end

end

