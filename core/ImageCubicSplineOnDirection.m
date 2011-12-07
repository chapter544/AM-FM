function y = ImageCubicSplineOnDirection(x, direction)
% Cubic Spline interpolation of a 2-D image
% function y = CubicSpline(x,K)
%
% Inputs:
%   @ x     : 2-D array input image
%   @ direction     : on row or on column (1 for row, others for column)
%
% Outputs: 
%   @ y     : cubic spline interpolation of input image x
%
% This algorithm is taken from M. Unser paper. The function will 
% perform a seperable operations (CubicSpline on row and then Cubic
% Spline on column).
%
% NTC - 3/2/2006

% Obtain row and column of image
m_row = size(x,1);
m_col = size(x,2);

% Allocate memory for output image
y = zeros(m_row,m_col);

if( direction == 1 ) % Perform CubicSpline on row 
	for m=1:1:m_row,
		y(m,:) = symExpFilt(x(m,:), -2+sqrt(3), 6.0);
	end
else % Perform CubicSpline on columns
	for n=1:1:m_col,
		y(:,n) = (symExpFilt( x(:,n)', -2+sqrt(3), 6.0) )';
	end
end
end % end ImageCubicSplineOnDirection



	
function y = symExpFilt( inData, pole, gain )
	len = length(inData);
	if( len < 2 )
		error('Insufficient data');
	end

	% Allocate memory for cPlus and cMinus array
	cPlus = zeros(1,len);
	cMinus = zeros(1,len);

	% Initialize coefficient
	cPlus(1) = inData(1);
	poleRow = pole;
	for m=2:len,
		cPlus(1) = cPlus(1) + inData(m) * poleRow;
		poleRow = poleRow * pole;
	end
	for m=len-1:-1:2,
		cPlus(1) = cPlus(1) + inData(m) * poleRow;
		poleRow = poleRow * pole;
	end
	cPlus(1) = cPlus(1) / (1.0 - poleRow);

	% Causal recursion
	for m=2:len,
		cPlus(m) = inData(m) + pole * cPlus(m-1);
	end

	% Mirror Symmetry
	cMinus(len) = cPlus(len) + pole*cPlus(len-1);
	cMinus(len) = cMinus(len) * (-gain * pole) / (1.0 - pole*pole);

	% Anti-causal recursion
	for m=len-1:-1:1,
		cMinus(m) = pole * ( cMinus(m+1) - gain*cPlus(m) );
	end

	y = cMinus;
end % end symExpFilt
