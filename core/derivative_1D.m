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
