function [t, y, meta] = load_ACSDSG2014(mode)
%load_ACSDSG2014 Load the digitized data from DSG and ZT's 2014 Resource loading
%ACS paper
% Call: 
%  Three Modes: 'RNAdeg', 'MGApt', and 'deGFP'
% output: 
% t: time points, in minutes
% y: # time points X # experimental conditions. data in nM
% meta: metadata for the experimental conditions. 
% the rawdata file should be in the same directory as this file
% DATA description: 
% !TODO fill this
fp = mfilename('fullpath');
slashes = regexp(fp, '/');
filedir = fp(1:slashes(end)-1);
run([filedir '/rawdata_ACSDSG2014.m'])
switch mode
    case 'RNAdeg'
        initialConc = [37.5, 75, 150, 200, 600, 700, 800, 900, 1000];
        % interpolate 0 to 120 , in 6 min intervals
        t = 60*(0:6:120)';
        y = zeros(length(t), length(rnadeg));
        for i = 1:length(rnadeg)
            y(:,i) = interp1(60*rnadeg{i}(:,1), rnadeg{i}(:,2), t, 'spline', 'extrap');
            % rescale the values so that they start exactly at the dosed
            % numbers
            y(:,i) = initialConc(i)*y(:,i)/y(1, i);
            
        end
        
        meta = {'t is in minutes'; 
            'initial RNA conc were 37.5, 75, 150, 200, 600, 700, 800, 900, 1000'};        
        
        
    case 'MGapt'
        t = 60*(0:8:800)';
        y = zeros(length(t), length(mg));
        for i = 1:length(mg)
            y(:,i) = interp1(60*mg{i}(:,1), mg{i}(:,2), t, 'spline', 'extrap');
        end
        meta = {'t is in minutes'; 
            'DNA conc for the mgAptamer were 0.5, 2, 5, 20'};
    case 'deGFP'
                t = 60*(0:8:800)';
        y = zeros(length(t), length(gfp));
        for i = 1:length(gfp)
            y(:,i) = interp1(60*gfp{i}(:,1), gfp{i}(:,2), t, 'spline', 'extrap');
        end
        meta = {'t is in minutes'; 
            'DNA conc for the GFP were 0.2, 0.5, 1, 2, 5, 20'};
        
    case 'deGFP_deriv'
        t = 60*(0:8:800)';
        y = zeros(length(t), length(gfp));
        dgfp = cell(length(gfp),1);
        
        
        for i = 1:length(gfp)
            tv_temp = 60*gfp{i}(1:end-1,1);
            dgfp_temp = diff(gfp{i}(:,2))./diff(60*gfp{i}(:,1));
            dgfp{i} = [tv_temp dgfp_temp];
            y(:,i) = interp1(dgfp{i}(:,1), dgfp{i}(:,2), t,...
                'spline', 'extrap');
        end
        meta = {'t is in minutes'; 
            'DNA conc for the GFP were 0.2, 0.5, 1, 2, 5, 20'};    
    end
end

