clear all
close all
clc

rootName = 'lena';
M = 256;
N = 256;


im_org = file2image2('float', M, N, sprintf('Images/%sR.float', rootName));

noise = 0.0 .* randn(M,N);

%im_org = imread('adelson.jpg');
%im_org = rgb2gray(im_org);
%im_org = double(im_org);
%[M,N] = size(im_org);



im_dc = mean(mean(im_org));
im = im_org - im_dc;
im = im + noise;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ANALYSIS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
numLevels = 5;
numOrien = 8;
%[bandFilterDft, LoDft] = GenerateSteerablePyramid(M, N, numLevels, numOrien);
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
R = cell(numLevels,numOrien);
T = cell(numLevels,numOrien);
Pls = cell(numLevels,numOrien);
unwrapOptions = [300 0 5 5 1 1];
for levidx=1:numLevels,
	for orienidx=1:numOrien,
		[A{levidx,orienidx}, ...
			U{levidx,orienidx}, V{levidx,orienidx}, P{levidx,orienidx}] = ...
			phaseUnwrap(chanResponseR{levidx,orienidx}, ... 
							chanResponseI{levidx,orienidx}, unwrapOptions);

		R{levidx,orienidx} = sqrt(U{levidx,orienidx}.^2 + V{levidx,orienidx}.^2);		
	end
end
