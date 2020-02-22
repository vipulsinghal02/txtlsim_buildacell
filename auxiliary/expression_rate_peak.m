function wholeExpFile = expression_rate_peak(wholeExpFile)

if(isfield(wholeExpFile,'diffReporterCurve'))
    tmp = max(wholeExpFile.diffReporterCurve);
else
    error('no derivative found, run calculateDerivaties first!');
end

wholeExpFile.ratePeak = tmp;

end

