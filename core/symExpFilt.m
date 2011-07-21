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
