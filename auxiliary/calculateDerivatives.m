function [whole_exp_file] = calculateDerivatives(whole_exp_file,varargin)


if isfield(whole_exp_file,'FileName')
    fileName = whole_exp_file.FileName;
else
    fileName = '';
end

if isfield(whole_exp_file,'Data')
    % original whole_exp_file
    t_vec = whole_exp_file.t_vec(1:whole_exp_file.numOfValidReads1)/60;
else
    % merged exp_file
    t_vec = whole_exp_file.t_vec/60;
end


if isfield(whole_exp_file,'MgCurveFit')
   NumOfFits = size(whole_exp_file.MgCurveFit,2);
   for k = 1:NumOfFits % number of Curves
     diffR(:,k,1) = differentiate(whole_exp_file.MgCurveFit(k).ff,whole_exp_file.t_vec/60);
   end
end

if isfield(whole_exp_file,'ReporterCurveFit')
    NumOfFits = size(whole_exp_file.ReporterCurveFit);
    for l = 1:NumOfFits(2) % number of Channels
        for k = 1:NumOfFits(1) % number of Curves
            diffR(:,k,l+1) = differentiate(whole_exp_file.ReporterCurveFit(k,l).ff,whole_exp_file.t_vec/60);
        end
    end
end

whole_exp_file.diffReporterCurve = diffR;


end