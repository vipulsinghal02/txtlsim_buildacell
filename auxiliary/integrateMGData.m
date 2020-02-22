function [whole_exp_file] = integrateMGData(whole_exp_file,varargin)

if nargin > 1
    fittypeName = varargin{1};
    f2 = fittype(fittypeName);
    opts = fitoptions( f2 );
else
    
    f2 = fittype( 'smoothingspline' );
    opts = fitoptions( f2 );
    opts.SmoothingParam = 2.5e-05;
end

% search for mg channel
channel_index = find(cellfun(@(x) strcmpi(x,'MG'),whole_exp_file.channels(:,1)) > 0);


if isfield(whole_exp_file,'FileName')
    fileName = whole_exp_file.FileName;
else
    fileName = '';
end

if isfield(whole_exp_file,'Data')
    % original whole_exp_file
    dataMtx = whole_exp_file.noBg(:,:,channel_index);
else
    % merged exp_file
    dataMtx = whole_exp_file.noBg_mean(:,:,channel_index);
end

t_vec = whole_exp_file.t_vec(1:size(dataMtx,1))/60;


disp(sprintf('calculating the int(MG) for individual wells in %s...',fileName));

tic
interval = 5;
intPoints = round(size(t_vec,1)/interval);
for k = 1:intPoints
    if interval*k > size(t_vec,1);
        idx = 1:size(t_vec,1);
    else 
        idx = 1:interval*k;
    end
    intMGData(k,:) = sum(dataMtx(idx,:));
end
toc
whole_exp_file.intMGData = intMGData;



end