function [A, U, V, P, Pls] = AMFM_TransformNoDC(im, numLevels, numOrien, varargin)

if( nargin < 3 )
	error('Not enough input argument: AMFM_Transform(im, numLevels, numOrien)');
end

% default values
unwrapOptions = [1 0 5 5 0.5];
if( length(varargin) == 1 )
	unwrapOptions = varargin{1};
end


% eliminate the DC portion
im = im - mean(im(:));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ANALYSIS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[bandFilterDft] = GenerateSteerablePyramidNoDC(size(im,1), size(im,2), numLevels, numOrien);

chanResponseR = cell(numLevels, numOrien);
chanResponseI = cell(numLevels, numOrien);

% Perform filtering through steerable pyramid
imDft = fftshift(fft2(im));
out = zeros(size(im));
for levidx = 1:numLevels,
	for orienidx = 1:numOrien,
		% Get the real part
		outDft = imDft .* bandFilterDft{levidx,orienidx};		
		chanResponseR{levidx,orienidx} = real(ifft2(ifftshift(outDft)));

		% Construct the imaginary part with the partial Hilbert transform
		H = GenerateDirectionalHilbertFilter(size(im,1), size(im,2), orienidx, numOrien);
		outImagDft = -sqrt(-1) .* outDft .* H; 
		chanResponseI{levidx,orienidx} = real(ifft2(ifftshift(outImagDft)));
	end
end
%Lo = real(ifft2(ifftshift(LoDft.*imDft)));


%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% DEMODUALTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
A = cell(numLevels,numOrien);
U = cell(numLevels,numOrien);
V = cell(numLevels,numOrien);
P = cell(numLevels,numOrien);
Pls = cell(numLevels,numOrien);
for levidx=1:numLevels,
	for orienidx=1:numOrien,
		[A{levidx,orienidx}, ...
			U{levidx,orienidx}, V{levidx,orienidx}, ...
			P{levidx,orienidx}, Pls{levidx,orienidx}] = ...
			phaseUnwrap(chanResponseR{levidx,orienidx}, ... 
							chanResponseI{levidx,orienidx}, unwrapOptions);
	end
end


