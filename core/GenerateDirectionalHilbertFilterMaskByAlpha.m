% Generate Discrete Directional Partial Hilbert transform kernel spectrum mask
% function [H] = GenerateDirectionalHilbertFilterMaskByAlpha(row, col, alpha)
%
% 	inputs:
%			row: vertical dimension
%			col: vertical dimension
%			alpha: rotation angle interm of radiant
%
%	return:
%			H: multiplication mask (either +1 or -1, or 0)
%
function [H] = GenerateDirectionalHilbertFilterMaskByAlpha(row, col, alpha)

H = zeros(row,col);
halfRow = ceil( (row-1) / 2 );
halfCol = ceil( (col-1) / 2 );


%[mm,nn] = meshgrid(-halfRow:1:halfRow-1,-halfCol:1:halfCol-1);
[mm,nn] = meshgrid(-halfCol:1:halfCol-1, -halfRow:1:halfRow-1);
H = sign( mm*cos(alpha) + nn*sin(alpha) );

% fix up the DC part
H(halfRow+1,halfCol+1) = 0.0;

% forces 0 at (-N/2,-N/2), (0,-N/2), and (-N/2,0)
H(halfRow+1,1) = 0.0;
H(1,1) = 0.0;
H(1,halfCol+1) = 0.0;

% fix up the correspondent on the N/2 lines
H(2:halfRow,1) = -1;
H(halfRow+2:end,1) = 1;
H(1,2:halfCol) = -1;
H(1,halfCol+2:end) = 1;



% replaced by the meshgrid stuff
%for m=2:row,
%	for n=2:col,
%		gain = (m-halfRow-1) * sin(alpha) + (n-halfCol-1) * cos(alpha);
%
%		if( gain > 0 )
%			H(m,n) = 1.0;
%		elseif (gain < 0) 
%			H(m,n) = -1.0;
%		else
%			if( (m-halfRow-1) > 0 )
%				H(m,n) = 1.0;
%			else
%				H(m,n) = -1.0;
%			end
%		end
%	end
%end


