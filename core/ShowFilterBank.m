% Generate "pseudo" frequency response of the steerable pyramid
%
% function [G] = ShowFilterBank(M,N, numLevel, numOrien)
%
%	G: G is the "pseudo" DFT of the steerable pyramid
%
function [G] =  ShowFilterBank(M,N, numLevel, numOrien)
[bandFilterDft] = GenerateSteerablePyramidNoDC(M,N, numLevel, numOrien);

G = zeros(M,N);
out = zeros(M,N);
for idx=1:numLevel,
	for idx2=1:numOrien,
		temp = abs(bandFilterDft{idx,idx2});
		G = max( G , temp );
	end
end
%image2file(out, 'float', 'steerablePyramid.float', 0);
