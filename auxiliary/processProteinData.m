function wholeExpFile = processProteinData(wholeExpFile)



% search for the protein reporter channels
channel_idx = find(cellfun(@(x) ~strcmpi(x,'MG'),wholeExpFile.channels(:,1)) > 0);

if isfield(wholeExpFile,'FileName')
    fileName = wholeExpFile.FileName;
else
    fileName = '';
end

if isfield(wholeExpFile,'Data')
    % original wholeExpFile
    dataMtx = wholeExpFile.noBg;
else
    % merged exp_file
    dataMtx = wholeExpFile.noBg_mean;
end

t_vec = wholeExpFile.t_vec(1:size(dataMtx,1))/60;

disp(sprintf('calculating the protein reporter curves fit for individual wells in %s...',fileName));
tic
for j=1:size(channel_idx,1)
    for k = 1:size(dataMtx,2);
        f2 = fittype( 'smoothingspline' );
        opts = fitoptions( f2 );
        if dataMtx(end,k,channel_idx(j)) < 3000
            opts.SmoothingParam = 1e-4;
        else
            opts.SmoothingParam = 0.0003039;
        end
        
        
        [fitinfo2(k,j).ff gof2(k,j).gof] = fit(t_vec,dataMtx(:,k,channel_idx(j)),f2,opts);
        wholeExpFile.ReporterCurve(:,k,j) = fitinfo2(k).ff(wholeExpFile.t_vec/60);
    end
end
toc
wholeExpFile.ReporterCurveFit = fitinfo2;
wholeExpFile.ReporterCurveFitGOF = gof2;

end