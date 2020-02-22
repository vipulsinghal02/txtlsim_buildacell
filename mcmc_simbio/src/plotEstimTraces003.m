function [ output_args ] = plotEstimTraces003(m,exportedMdlObj,tspan, ...
    simulatedDataMatrix, dosedInitVals,...
    measuredSpecies, varargin)
%plotEstimTraces003 Plot estimated trace mean and standard deviation for
%txtl data
%   

p = inputParser;
p.addParameter('paramID',[10 12 13],@isnumeric);
p.addParameter('titlestr',[] ,@ischar)
p.parse(varargin{:});
p=p.Results;

colIdx = p.paramID;
if ndims(m) == 3
    m = m(:,:)';
end
    
I = sampleIntersections(m, colIdx, 'mass', 0.4);

numPostPlots = min([250, length(I)]); %plot numPostPlots samples...
if numPostPlots== 0
    error('No parameter points meet specified search criteria')
elseif numPostPlots<30
    warning(['Only using ' num2str(numPostPlots) 'parameter points. Expand Search parameters'])
end

r = ceil(rand(numPostPlots,1)*size(I,1));
% I = (1:size(msimple,1))';
idx = I(r);
conc = simulatedDataMatrix;
nICs = size(dosedInitVals,1);
if length(measuredSpecies)==2
[measuredNames{1}, measuredNames{2}] = deal(measuredSpecies.objectName);
elseif length(measuredSpecies)==1
    [measuredNames{1}] = deal(measuredSpecies.objectName);
end



% Compute sample trajectories, max values for axis limits, means, standard deviations 
meanGFP = zeros(length(tspan), nICs);
meanRNA = zeros(length(tspan), nICs);
stdevGFP = zeros(length(tspan), nICs);
stdevRNA = zeros(length(tspan), nICs);
GFPsampleTraj = zeros(length(tspan), nICs, numPostPlots);
RNAsampleTraj = zeros(length(tspan), nICs, numPostPlots);

idxnotused = [];
kknotused = [];
for i = 1:nICs
    for kk=1:numPostPlots
        try
        sd = simulate(exportedMdlObj, [exp(m(idx(kk),:)'); dosedInitVals(i,:)']);
        sd = resample(sd, tspan);
        spSD = selectbyname(sd, measuredNames); 
        [RNAsampleTraj(:,i, kk), GFPsampleTraj(:,i,kk)] = deal(spSD.Data(:,1), spSD.Data(:,2));
        catch
            idx(kk)
            idxnotused = [idxnotused; idx(kk)];
            kknotused = [kknotused;kk];
            
        end
    end
    % compute mean and std
    
    meanGFP(:,i)=mean(GFPsampleTraj(:,i,setdiff(1:numPostPlots, kknotused)),3);
    meanRNA(:,i)=mean(RNAsampleTraj(:,i,setdiff(1:numPostPlots, kknotused)),3);
    stdevGFP(:,i)=std(GFPsampleTraj(:,i,setdiff(1:numPostPlots, kknotused)),0,3);
    stdevRNA(:,i)=std(RNAsampleTraj(:,i,setdiff(1:numPostPlots, kknotused)),0,3);
end

%%

%
cc = colorschemes;
figure
ss = get(0, 'screensize');
set(gcf, 'Position', [50 100 ss(3)/1.1 ss(4)/1.3]);
%
    maxRNA = max(max(max(meanRNA+stdevRNA)), max(conc(:,1)));
    maxGFP = max(max(max(meanGFP+stdevGFP)), max(conc(:,2)));
    
for i = 1:nICs
    
    % RNA
    subplot(nICs, 2,2*(i-1)+1);
    [h(1), ptch(1)] = boundedline(tspan/3600, meanRNA(:,i), stdevRNA(:,i));
    set(ptch(1), 'FaceColor', cc{2,7}(1,:), 'FaceAlpha', 0.5);
    set(h(1), 'Color', cc{2,7}(2,:), 'LineStyle', '--');
    hold on 
    set(h(1), 'LineWidth', 2)
    h(2)=plot(tspan/3600,conc((i-1)*161 + 1 : (i)*161,1),'color',cc{2,7}(3,:) ,'linewidth',2);
    hold on 
    set(gca, 'Ylim', [0, round(maxRNA+5, -1)])
    xlabel('time/h')
    ylabel('conc/nM')
    title(sprintf('%s conc, DNA=%0.2g nM', 'RNA', dosedInitVals(i,:)'))
    
    % GFP
    subplot(nICs, 2,2*(i-1)+2);
    [h(3), ptch(2)] = boundedline(tspan/3600, meanGFP(:,i), stdevGFP(:,i));
    set(ptch(2), 'FaceColor', cc{2,9}(3,:), 'FaceAlpha', 0.5);
    set(h(3), 'Color', cc{2,9}(2,:), 'LineStyle', '--');
    hold on 
    set(h(3), 'LineWidth', 2)
    h(4)=plot(tspan/3600,conc((i-1)*161 + 1 : (i)*161,2),'color',cc{2,9}(1,:) ,'linewidth',2);
    hold on 
    set(gca, 'Ylim', [0, round(maxGFP+10, -1)])
    xlabel('time/h')
    ylabel('conc/nM')
    title(sprintf('%s conc, DNA=%0.2g nM', 'GFP', dosedInitVals(i,:)'))
    
    
    ax = gca;
end
if ~isempty(p.titlestr)
suptitle(p.titlestr)
end
 legend(ax, [h(1), h(2), h(3), h(4)],{'RNA sample mean','RNA Truth',...
     'GFP sample mean',...
      'GFP Truth'})

  %% Plot the true traces, estimated means and shade the standard deviations. 
  figure
ss = get(0, 'screensize');
set(gcf, 'Position', [50 100 ss(3)/1.1 ss(4)/1.3]);
%
    maxRNA = max(max(max(meanRNA+stdevRNA)), max(conc(:,1)));
    maxGFP = max(max(max(meanGFP+stdevGFP)), max(conc(:,2)));
    
for i = 1:nICs
    
    % RNA
    subplot(nICs, 2,2*(i-1)+1);
    [h(1), ptch(1)] = boundedline(tspan/3600, meanRNA(:,i), stdevRNA(:,i));
    set(ptch(1), 'FaceColor', cc{2,7}(1,:), 'FaceAlpha', 0.5);
    set(h(1), 'Color', cc{2,7}(2,:), 'LineStyle', '--');
    hold on 
    set(h(1), 'LineWidth', 2)
    h(2)=plot(tspan/3600,conc((i-1)*161 + 1 : (i)*161,1),'color',cc{2,7}(3,:) ,'linewidth',2);
    hold on 
    set(gca, 'Ylim', [0, round(maxRNA+5, -1)])
    xlabel('time/h')
    ylabel('conc/nM')
    title(sprintf('%s conc, DNA=%0.2g nM', 'RNA', dosedInitVals(i,:)'))
    
    % GFP
    subplot(nICs, 2,2*(i-1)+2);
    [h(3), ptch(2)] = boundedline(tspan/3600, meanGFP(:,i), stdevGFP(:,i));
    set(ptch(2), 'FaceColor', cc{2,9}(3,:), 'FaceAlpha', 0.5);
    set(h(3), 'Color', cc{2,9}(2,:), 'LineStyle', '--');
    hold on 
    set(h(3), 'LineWidth', 2)
    h(4)=plot(tspan/3600,conc((i-1)*161 + 1 : (i)*161,2),'color',cc{2,9}(1,:) ,'linewidth',2);
    hold on 
    set(gca, 'Ylim', [0, round(maxGFP+10, -1)])
    xlabel('time/h')
    ylabel('conc/nM')
    title(sprintf('%s conc, DNA=%0.2g nM', 'GFP', dosedInitVals(i,:)'))
    
    
    ax = gca;
end

if ~isempty(p.titlestr)
suptitle(p.titlestr)
end
 legend(ax, [h(1), h(2), h(3), h(4)],{'RNA sample mean','RNA Truth',...
     'GFP sample mean',...
      'GFP Truth'})

end

