function [dx, dy] = SmoothCubicSplineGradientApproximation(y, lamda)
% Smooth Cubic Spline Gradient Approximation
% function [dx, dy] = SmoothCubicSplineGradientApproximation(y, lamda)
%
% Inputs:
%   @ y     : 2-D array input image
%   @ lamda : variance 
%
% Outputs: 
%   @ dx     : Horizontal gradient
%   @ dy     : Vertical gradient
%
% This algorithm is taken from M. Unser paper, it is the appximation
% of the smooth spline (page 843, equation 4.9)
%%


smoothVariance = sqrt(2*lamda);

[M,N] = size(y);
h = fspecial('gaussian', [5 5], smoothVariance)
y_filt = imfilter(y, h, 'replicate');

% Allocate memory for output image
dx = zeros(size(y));
dy = zeros(size(y));

% Horizontal gradient
for m=1:M
	dx(m,:) = derivative_1D(y_filt(m,:));
end

% Vertical gradient
for n=1:N
	dy(:,n) = derivative_1D(y_filt(:,n)')';
end

