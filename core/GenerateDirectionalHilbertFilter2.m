function [H] = GenerateDirectionalHilbertFilter2(row, col, orienidx, numOrien)
out = zeros(row,col);
halfRow = row / 2;
halfCol = row / 2;
alpha = (pi/(2*numOrien)) + (orienidx - 1) / numOrien * pi - pi/2;

for m=2:row,
	for n=2:col,
		gain = (m-halfRow-1) * sin(alpha) + (n-halfCol-1) * cos(alpha);

		if( gain > 0 )
			out(m,n) = 1.0;
		elseif (gain < 0) 
			out(m,n) = -1.0;
		else
			if( (m-halfRow-1) > 0 )
				out(m,n) = 1.0;
			else
				out(m,n) = -1.0;
			end
		end
	end
end

% fix up the DC part
out(halfRow+1,halfCol+1) = 0.0;

% fix up the correspondent
out(:,1) = -out(:,halfCol+1);
out(1,:) = out(halfRow+1,:);

H = out;
