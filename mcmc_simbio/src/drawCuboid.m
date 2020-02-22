function [outputArg1,outputArg2] = drawCuboid(inputArg1,inputArg2)
% Need to work on this .Useful in script analysis_vnprl_F2.m
%   Detailed explanation goes here

xc=1; yc=1; zc=1;    % coordinated of the center
L=10;                 % cube size (length of an edge)
alpha=0.8;           % transparency (max=1=opaque)

X = [0 0 0 0 0 1; 1 0 1 1 1 1; 1 0 1 1 1 1; 0 0 0 0 0 1];
Y = [0 0 0 0 1 0; 0 1 0 0 1 1; 0 1 1 1 1 1; 0 0 1 1 1 0];
Z = [0 0 1 0 0 0; 0 0 1 0 0 0; 1 1 1 0 1 1; 1 1 1 0 1 1];

C='blue';                  % unicolor

X = L*(X-0.5) + xc;
Y = L*(Y-0.5) + yc;
Z = L*(Z-0.5) + zc; 

fill3(X,Y,Z,C,'FaceAlpha',alpha);    % draw cube
axis equal


end

