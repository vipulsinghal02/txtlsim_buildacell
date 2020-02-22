function [fig, ax] = mcmc_3D(mcut, axL, titl)
	% mcmc_3D takes an array of dimensions nPoints x 3 and generates a
	% scatterplot. 
	%
	% mcut 	is a nPoints x 3 matrix of points to scatterplot. 
	% axL 	must be a 1 x 3 cell array of the axis label strings 
	% 		corresponding to the 3 columns of mcus. 
	% titl 	is a title string. 
	% 
fig = figure;


XX = mcut(:, 1);
YY = mcut(:, 2);
ZZ = mcut(:, 3);
scatter3(XX,YY,ZZ)
xlabel(axL{1}, 'FontSize', 20)
ylabel(axL{2}, 'FontSize', 20)
zlabel(axL{3}, 'FontSize', 20)
title(titl, 'FontSize', 20)
ax = gca;
end
