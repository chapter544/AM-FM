%Perform AM-FM transform
% function [A, U, V, P, Pls, Resi] = 
%	AMFM_Transform(im, numLevels, numOrien, varargin)
%paramters:
%	im: 2-D input real image array
%	numLevels: number of levels for decomposition
%	numOrien: number of orientations for decomposition
%
%	Optional argument:
%	params	: a vector containg 5 elements [noDC PHASE_CONST smoothing L sigma]
%		noDC: 1 for withDCRing, 0 for without DC ring
%		phaseUnwrapping: 1 for perform phase unwrapping, 0 for NO
%		PHASE_CONSTANT: set to 300 for smoothness of the unwrapped phase
%		smoothing = 1 for Y, 0 for N)
%		L = width of smoothing window
%		sigma = bandwidth of smoothing filter (variance of gaussian)
%		
%		default: [1 1 300 1 3  0.5]
%	
%return:
%	A: AM function, in cell structure, i.e., A{1,1}, A{1,2}, A{3,1} ...
%	U: FM function in horizontal direction
%	V: FM function in vertical direction
%	P: unwrapped phase with congruency enforced
%	Pls: unwrapped phase directly from the 2-D least square phase 
%			unwrapping algorithm 
%	Resi: residual of the decomposition wrt the original input image
%
%	All return values are in structure types,
%		i.e., A{1,1}, A{1,2}, A{3,1} ...
function [A, U, V, P, Pls, Resi] = ...
	AMFM_Transform(im_org, numLevels, numOrien, varargin)

if( nargin < 3 )
	error('Not enough input argument: AMFM_Transform(im, numLevels, numOrien)');
end

if(size(im_org,3) > 1)
	error('Can only perform on 2D image');
end

if(numOrien < 1)
	error('numOrien has to be integer >= 1');
end
if(numLevels < 1)
	error('numLevels has to be integer >= 1');
end


% default values
trans_opts = [0 1 300 1 3 0.5];
if( length(varargin) == 1 )
	trans_opts = varargin{1};
else
	fprintf('Performing AM-FM transform using default parameters:\n');
	fprintf('-Exclude the low-frequency disk: %d\n',trans_opts(1));
	fprintf('-Demodulation parameters: \n');
	if(trans_opts(2) == 1) % phase unwrapping: YES, NO
		fprintf('\tPhase Unwrapping: YES\n');
		fprintf('\tPHASE_CONST: %0.2f\n',trans_opts(3));
	else
		fprintf('\tPhase Unwrapping: NO');
	end
	fprintf('\tsmoothing: %d\n',trans_opts(4));
	if(trans_opts(3) == 1) 
		fprintf('\tGaussian kernel\n');
		fprintf('\twindow: %dx%d\n',trans_opts(5), trans_opts(5));
		fprintf('\tsigma: %0.2f\n',trans_opts(6));
	end
end

% eliminate the DC portion
imDC = mean(im_org(:));
im = im_org - imDC;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ANALYSIS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
imDft = fftshift(fft2(im));
if(trans_opts(1) == 1) % with DC Ring as residual
	[bandFilterDft, LoDft] = GenerateSteerablePyramid(...
		size(im,1), size(im,2), numLevels, numOrien);

	Resi = real(ifft2(ifftshift(LoDft.*imDft))) + imDC;
else
	[bandFilterDft] = GenerateSteerablePyramidNoDC(...
		size(im,1), size(im,2), numLevels, numOrien);

	Resi = imDC;
end

% Perform filtering through steerable pyramid
chanResponseR = cell(numLevels, numOrien);
chanResponseI = cell(numLevels, numOrien);
out = zeros(size(im));
for levidx = 1:numLevels,
	for orienidx = 1:numOrien,
		% Get the real part
		outDft = imDft .* bandFilterDft{levidx,orienidx};		
		chanResponseR{levidx,orienidx} = real(ifft2(ifftshift(outDft)));

		% Construct the imaginary part with the partial Hilbert transform
		%H = DirectionalHilbertFilterMaskByOrientation(...
		%	size(im,1), size(im,2), orienidx, numOrien);

		alpha = (pi/(2*numOrien)) + (orienidx - 1) / numOrien * pi - pi/2;
		H = GenerateDirectionalHilbertFilterMaskByAlpha(...
			size(im,1), size(im,2), alpha);
		outImagDft = -sqrt(-1) .* outDft .* H; 
		chanResponseI{levidx,orienidx} = real(ifft2(ifftshift(outImagDft)));
		
%		% create by alpha mask with adjusted points 
%		alpha = (pi/(2*numOrien)) + (orienidx - 1) / numOrien * pi - pi/2;
%		H2 = GenerateDirectionalHilbertFilterMaskByAlpha(...
%			size(im,1), size(im,2), alpha);
%		outImagDft2 = -sqrt(-1) .* outDft .* H2; 
%		blah = real(ifft2(ifftshift(outImagDft2)));
%
%
%		a = sqrt(chanResponseI{levidx,orienidx}.^2 + chanResponseR{levidx,orienidx}.^2);
%		a2 = sqrt(blah.^2 + chanResponseR{levidx,orienidx}.^2);
%
%		diff = blah - chanResponseI{levidx,orienidx};
%		diffA = a - a2;
%
%		subplot(3,2,1)
%		imagesc(blah);axis('image'); axis off; colormap(gray);
%		subplot(3,2,2)
%		imagesc(chanResponseI{levidx,orienidx});axis('image'); axis off; colormap(gray);
%		subplot(3,2,3)
%		imagesc(a2);axis('image'); axis off; colormap(gray);
%		title('H adjusted');
%		subplot(3,2,4)
%		imagesc(a);axis('image'); axis off; colormap(gray);
%		title('H org');
%		subplot(3,2,5)
%		imagesc(diff);axis('image'); axis off; colormap(gray);
%		title(sprintf('min=%f, max=%f', min(min(diff)), max(max(diff))));
%		subplot(3,2,6)
%		imagesc(diffA);axis('image'); axis off; colormap(gray);
%		title(sprintf('min=%f, max=%f', min(min(diffA)), max(max(diffA))));
%
%		pause;
	end
end
%---------------------------------------------------------------------------%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% DEMODUALTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
				chanResponseI{levidx,orienidx}, trans_opts(2:end));
	end
end
%------------------------------------------------------------------------%
