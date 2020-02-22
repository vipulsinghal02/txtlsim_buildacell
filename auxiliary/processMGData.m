function [whole_exp_file] = processMGData(whole_exp_file,varargin)

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


disp(sprintf('calculating the MG curve fit for individual wells in %s...',fileName));
tic
for k = 1:size(dataMtx,2);
    [fitinfo2(k).ff gof2(k).gof] = fit(t_vec,dataMtx(:,k),f2,opts);
    whole_exp_file.MgCurve(:,k) = fitinfo2(k).ff(whole_exp_file.t_vec/60);
end
toc
whole_exp_file.MgCurveFit = fitinfo2;
whole_exp_file.MgCurveFitGOF = gof2;




end