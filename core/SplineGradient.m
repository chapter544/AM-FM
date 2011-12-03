function [dx, dy] = SplineGradient(y)
% Cubic Spline Gradient
% function [dx, dy] = SplineGradient(y)
%
% Inputs:
%   @ y     : 2-D array input image
%
% Outputs: 
%   @ dx     : Horizontal gradient
%   @ dy     : Vertical gradient
%
% This algorithm is taken from M. Unser paper. The function will 
% perform a seperable operations (CubicSpline on row and then Cubic
% Spline on column).
%
% NTC - 3/2/2006

% Obtain row and column of image
row = size(y,1);
col = size(y,2);

% Allocate memory for output image
dx = zeros(size(y));
dy = zeros(size(y));

im_H = ImageCubicSplineOnDirection(y, 1);
im_V = ImageCubicSplineOnDirection(y, 2);

% Horizontal gradient
for m=1:row,
	dx(m,:) = derivative_1D(im_H(m,:));
end

% Vertical gradient
for n=1:col,
	dy(:,n) = derivative_1D(im_V(:,n)')';
end


% derivative_1D Calculate the discrete derivative with
% kernel [0.5 0 -0.5]
% 
% NTC 06/19/08
%
function y = derivative_1D(x)

% find the length of array
N = length(x);
y = zeros(1,N);

for m=2:N-1,
	y(m) = 0.5*(x(m+1) - x(m-1));
end

end %end derivative_1D
end %end SplineGradient
