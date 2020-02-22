function h = stdshade(xVector,meanData,stdData,varargin)
% usage: stdshade(xVector,meanData,stdData,alpha,acolor)
% the function plots the meanData vector/matrix, at which each column is an observation. 
% the corresponding stdData is shown as shading.
% - xVector contains the values of the x-axis
% - meanData column vector of data
% - stdData column vector of data used for shadowing, typically std of the
% data 
% - alpha defines transparency of the shading
% - acolor defines the used color (default is red)

assert(all(size(meanData) == size(stdData)))

if nargin == 3
    alpha = 0.3;
    acolor = 'r';
elseif nargin == 5
    alpha = varargin{1};
    acolor = varargin{2};
end

% more than one data set is given
numOfDataSets = size(meanData,2);
if  numOfDataSets > 1
    if size(alpha,1) ~= numOfDataSets
        alpha = repmat(alpha,numOfDataSets,1);
    end
    for k=1:numOfDataSets
        h(k) = stdshade(xVector,meanData(:,k),stdData(:,k),alpha(k),acolor{k});
    end
else

    fill([xVector' fliplr(xVector')],[meanData'+stdData' fliplr(meanData'-stdData')],acolor, 'FaceAlpha', alpha,'linestyle','none');
    
    if ishold==0
        check=true;
    else
        check=false;
    end
    
    hold on;
    h = plot(xVector,meanData,'Color',acolor,'linewidth',1.5); %% change color or linewidth to adjust mean line
    
    if check
        hold off;
    end
    
end

end



