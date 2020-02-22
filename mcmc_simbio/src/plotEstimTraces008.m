function plotEstimTraces008(m, genemodel, ICarray, tspan,...
    spconc1, spconc2, varargin)
%plotEstimTraces008 Sample from the posterior distributions, and generate
%time traces, means, standard deviations. 
% Optional Capability: if true parameters are provided, plot the true data too. 
% 
%
%
%
%
nSp = 3;
if ndims(m) == 3
    m = m(:,:)';
end

p = inputParser;
p.addParameter('mode','samplesonly',@ischar); % 'samplesonly', 'trueparams'
p.addParameter('sampling', 'all', @ischar) % 'all', 'MAPint, 'percentiles'
p.addParameter('paramID',[10 12 13],@isnumeric);
p.addParameter('titlestr',[] ,@ischar)
p.addParameter('ncurves',200 ,@isnumeric)
p.parse(varargin{:});
p=p.Results;

ncurves = p.ncurves; 

if strcmp(p.sampling, 'all')
    idx = ceil(rand(ncurves, 1)*size(m,1));
elseif strcmp(p.sampling, 'percentiles')
    %!TODO, percentiles
elseif strcmp(p.sampling, 'MAPint')
    %!TODO: MAP intervals, some mass around
end

nICs = size(ICarray,1);
% compute the mean and the axis limits

estimconc1 = zeros(length(tspan),nSp, nICs , ncurves);
estimconc2 = zeros(length(tspan),nSp, nICs , ncurves);
maxRNA_sd = 0;
maxGFP_sd = 0;
for i = 1:nICs
    for kk=1:ncurves
        sp0 = ICarray(i,:);
        [~,estimconc1(:,:, i, kk)] = genemodel(m(idx(kk),[1:4 7 8]), sp0, tspan);
        [~,estimconc2(:,:, i, kk)] = genemodel(m(idx(kk),[1:2 5:8]), sp0, tspan); 
    end
    % compute mean and std
    mean1 = mean(estimconc1, 4); % mean in dim 4
    mean2 = mean(estimconc2, 4);
    std1 = std(estimconc1, 0, 4); % SD in dim 4
    std2 = std(estimconc2, 0, 4); 
    maxRNA_sd = max([maxRNA_sd, max(mean1(:,2, i)+std1(:,2, i)), max(mean2(:,2,i)+std2(:,2,i))]);
    maxGFP_sd = max([maxGFP_sd, max(mean1(:,3,i)+std1(:,3,i)), max(mean2(:,3,i)+std2(:,3,i))]);
end

maxRNA_curves = max(max(max(max([estimconc1(:,2,:,:), estimconc2(:,2,:,:)]))));
maxGFP_curves = max(max(max(max([estimconc1(:,3,:,:), estimconc2(:,3,:,:)]))));

cc = colorschemes;
% figure
% ss = get(0, 'screensize');
% set(gcf, 'Position', [50 100 ss(3)/1.1 ss(4)/1.3]);
% % plot the mean and individual curves
% for i = 1:nICs
%     for kk=1:ncurves
%         subplot(nICs, 4,4*(i-1)+1);
%         h(1)=plot(tspan,estimconc1(:,2,i,kk),'color',[.6 .35 .3].^.2);
%         hold on
%         subplot(nICs, 4,4*(i-1)+2); 
%         h(2)=plot(tspan,estimconc1(:,3,i,kk),'color',[.2 .75 .2].^.2);
%         hold on
%         subplot(nICs, 4,4*(i-1)+3); 
%         h(3)=plot(tspan,estimconc2(:,2,i,kk),'color',[.6 .35 .3].^.2);
%         hold on
%         subplot(nICs, 4,4*(i-1)+4); 
%         h(4)=plot(tspan,estimconc2(:,3,i,kk),'color',[.2 .75 .2].^.2);
%         hold on        
%     end
%     subplot(nICs, 4,4*(i-1)+1);
%     h(5)=plot(tspan,spconc1(:,2,i),'r','linewidth',2);
%     h(9) = plot(tspan,mean1(:,2,i),'color',[.6 .35 .3],'linewidth',2, 'LineStyle',':'); 
%     title(sprintf('%s, E%d, DNA = %0.2g', 'RNA', 1, ICarray(i,1)))
%     set(gca, 'Ylim', [0, round(maxRNA_curves+5)])
%     hold on 
%     subplot(nICs, 4,4*(i-1)+2); 
%     h(6)=plot(tspan,spconc1(:,3,i),'g','linewidth',2);  
%     h(10) = plot(tspan,mean1(:,3,i),'color',[.2 .75 .2],'linewidth',2, 'LineStyle',':'); 
%     title(sprintf('%s, E%d, DNA = %0.2g', 'GFP', 1, ICarray(i,1)))
%     set(gca, 'Ylim', [0, round(maxGFP_curves+5)])
%     hold on 
%     subplot(nICs, 4,4*(i-1)+3);
%     h(7)=plot(tspan,spconc2(:,2,i),'r','linewidth',2);
%     h(11) = plot(tspan,mean2(:,2,i),'color',[.6 .35 .3],'linewidth',2, 'LineStyle',':');
%     title(sprintf('%s, E%d, DNA = %0.2g', 'RNA', 2, ICarray(i,1)))
%     set(gca, 'Ylim', [0, round(maxRNA_curves+5)])
%     hold on 
%     subplot(nICs, 4,4*(i-1)+4); 
%     h(8)=plot(tspan,spconc2(:,3,i),'g','linewidth',2); 
%     h(12) = plot(tspan,mean2(:,3,i),'color',[.2 .75 .2],'linewidth',2, 'LineStyle',':');
%     title(sprintf('%s, E%d, DNA = %0.2g', 'DNA', 2, ICarray(i,1)))
%     set(gca, 'Ylim', [0, round(maxGFP_curves+5)])
%     hold on 
% end
% 
% axis tight
% legend([h(1), h(2), h(5), h(6), h(9), h(10)],'RNA Samples','GFP Samples',...
%     'RNA True', 'GFP True', 'RNA mean', 'GFP mean')
% suplabel('Plots using parameter estimates'  ,'t');

% plot the mean and standard deviation
figure
ss = get(0, 'screensize');
set(gcf, 'Position', [50 100 ss(3)/1.8 ss(4)/1.8]);
for i = 1:nICs
    % compute mean and std
    
    subplot(nICs, 4,4*(i-1)+1);
    h(5)=plot(tspan,spconc1(:,2,i),'color',cc{2,7}(2,:) ,'linewidth',2);
    hold on 
    [h(1), ptch(1)] = boundedline(tspan, mean1(:,2,i), std1(:,2,i));
    set(ptch(1), 'FaceColor', cc{2,7}(1,:), 'FaceAlpha', 0.5);
    set(h(1), 'Color', cc{2,7}(2,:).^2, 'LineStyle', '--');
    hold on 
    set(h(1), 'LineWidth', 2)

    set(gca, 'Ylim', [0, round(maxRNA_sd+5)])
    title(sprintf('%s, E%d, DNA = %0.2g', 'RNA', 1, ICarray(i,1)), 'FontSize', 16)
    
    subplot(nICs, 4,4*(i-1)+2); 
    h(6)=plot(tspan,spconc1(:,3,i),'color',cc{2,9}(2,:),'linewidth',2);  
    title(sprintf('%s, E%d, DNA = %0.2g', 'GFP', 1, ICarray(i,1)), 'FontSize', 16)
    hold on
    [h(2), ptch(2)] = boundedline(tspan, mean1(:,3,i), std1(:,3,i));
    set(ptch(2), 'FaceColor', cc{2,9}(3,:), 'FaceAlpha', 0.5);
    set(h(2), 'Color', cc{2,9}(2,:).^2, 'LineStyle', '--');    
    hold on 
    set(h(2), 'LineWidth', 2)
 
    set(gca, 'Ylim', [0, round(maxGFP_sd+5)])
    
    subplot(nICs, 4,4*(i-1)+3);
    h(7)=plot(tspan,spconc2(:,2,i),'color',cc{2,7}(2,:) ,'linewidth',2);
    title(sprintf('%s, E%d, DNA = %0.2g', 'RNA', 2, ICarray(i,1)), 'FontSize', 16)
    hold on 
    [h(3),ptch(3)] = boundedline(tspan, mean2(:,2,i), std2(:,2,i));
    set(ptch(3), 'FaceColor', cc{2,7}(1,:), 'FaceAlpha', 0.5);
    set(h(3), 'Color', cc{2,7}(2,:).^2, 'LineStyle', '--');
    hold on 
    set(h(3), 'LineWidth', 2)

    set(gca, 'Ylim', [0, round(maxRNA_sd+5)])
   
    
    subplot(nICs, 4,4*(i-1)+4); 
    h(8)=plot(tspan,spconc2(:,3,i),'color',cc{2,9}(2,:),'linewidth',2); 
    title(sprintf('%s, E%d, DNA = %0.2g', 'GFP', 2, ICarray(i,1)), 'FontSize', 16)
    hold on 
    [h(4), ptch(4)] = boundedline(tspan, mean2(:,3,i), std2(:,3,i));
    set(ptch(4), 'FaceColor', cc{2,9}(3,:), 'FaceAlpha', 0.5);
    set(h(4), 'Color', cc{2,9}(2,:).^2, 'LineStyle', '--');     
    hold on 
    set(h(4), 'LineWidth', 2)

    set(gca, 'Ylim', [0, round(maxGFP_sd+5)])
    if i ==1
    ax = gca;
    end
end

%
%  axis tight
%   suplabel('Plots using parameter estimates'  ,'t');
  [ax1,h1] =suplabel('time, arbitrary units'  );
  [ax2,h2] =suplabel('concentration, arbitrary units'  ,'y');
 set(h1,'FontSize',20)
 set(h2,'FontSize',20)
 legend(ax, [h(3), h(4), h(7), h(8)],{'RNA mean','GFP mean',...
     'RNA True', 'GFP True'}, 'FontSize', 12, 'Location', 'NorthEast')



end
