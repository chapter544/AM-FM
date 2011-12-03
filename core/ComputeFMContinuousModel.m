% Compute the A, U, V without phase unwrapping using the continuous 
% demodulation algorithm.
% function [A U V] = computeFMContinuousModel(realImg, imagImg, varargin)
% 
% inputs:
%		realImg:	real input image
%		imagImg:	imaginary image
%
% outputs:
%		A:	AM function
%		U:	horizontal component of FM
%		V:	vertical component of FM
function [A, U, V] = ComputeFMContinuousModel(realImg, imagImg, varargin)

% Argument checking
if(nargin ~= 2)
	error('Require 2 input arguments: realImg and imagImg');
end

% check for image dimension argreement
if( size(realImg,1)  ~= size(imagImg,1) || ... 
	size(realImg,2)  ~= size(imagImg,2) )
	error('Real and imaginary image are of different size');
end

% default smoothing to 1 (always perform smoothing)
L = 0;
sigma = 0.5;
if(size(varargin,2) == 1)
	L = varargin{1};
elseif(size(varargin,2) == 2)
	L = varargin{1};
	sigma = varargin{2};
end

% specify the threshold value
LOW_THRESHOLD = 10^(-3);

% AM as the magnitude of the complex image
A = sqrt((realImg.*realImg + imagImg.*imagImg));
Asquare = A .* A;
Asquare = Asquare + LOW_THRESHOLD .* (Asquare < LOW_THRESHOLD); 


% Find gradient images
[realDev_H, realDev_V] = SplineGradient(realImg);
[imagDev_H, imagDev_V] = SplineGradient(imagImg);

% phase derivative (no unwrap phase needed).
U = (imagDev_H.*realImg - realDev_H.*imagImg) ./ Asquare;
V = (imagDev_V.*realImg - realDev_V.*imagImg) ./ Asquare;

% default gaussian kernel of size [3 3] and sigma = 0.5
if( L > 0 )
	h = fspecial('gaussian', sigma, [L L]);
	U = imfilter(U, h, 'replicate');
	V = imfilter(V, h, 'replicate');
end
