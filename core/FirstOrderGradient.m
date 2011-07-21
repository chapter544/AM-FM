function [dx, dy] = FirstOrderGradient(y)
row = size(y,1);
col = size(y,2);

% Allocate memory for output image
dx = zeros(size(y));
dy = zeros(size(y));

% Horizontal gradient
for m=1:row,
	for n=1:col-1,
		dx(m,n) = y(m,n+1) - y(m,n);
	end
end

% Vertical gradient
for n=1:col,
	for m=1:row-1,
		dy(m,n) = y(m+1,n) - y(m,n);
	end
end

