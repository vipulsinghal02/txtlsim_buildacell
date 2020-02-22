function [ fixedExpFile ] = replaceMissingData( sourceExpFile, destExpFile, sourceIndex, destIndex )

fixedExpFile = destExpFile;

% fields to fix:
% 
% Data: [281x23x3 double]
% noBg: [281x23x3 double]
% rate: {23x3 cell}
% rateFitInfo: {23x3 cell}
% endTimeIndex: {23x3 cell}
% endTime: {23x3 cell}
% MgCurve: [281x23 double]
% MgCurveFit: [1x23 struct]
% MgCurveFitGOF: [1x23 struct]

vectorLength = min(sourceExpFile.numOfValidReads,destExpFile.numOfValidReads);

fixedExpFile.Data(:,destIndex,:) = sourceExpFile.Data(:,sourceIndex,:);
fixedExpFile.noBg(:,destIndex,:) = sourceExpFile.noBg(:,sourceIndex,:);
fixedExpFile.rate(destIndex,:) = sourceExpFile.rate(sourceIndex,:);
fixedExpFile.rateFitInfo(destIndex,:) = sourceExpFile.rateFitInfo(sourceIndex,:);
fixedExpFile.endTimeIndex(destIndex,:) = sourceExpFile.endTimeIndex(sourceIndex,:);
fixedExpFile.endTime(destIndex,:) = sourceExpFile.endTime(sourceIndex,:);
fixedExpFile.MgCurve(:,destIndex) = sourceExpFile.MgCurve(:,sourceIndex);
fixedExpFile.MgCurveFit(destIndex) = sourceExpFile.MgCurveFit(sourceIndex);
fixedExpFile.MgCurveFitGOF(destIndex) = sourceExpFile.MgCurveFitGOF(sourceIndex);


end

