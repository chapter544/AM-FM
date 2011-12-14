%DiscreteDemodulation:  2-D  estimation of the AM and FM functions
%paramters:
%	realImg	: 2-D input real image array
%	imagImg	: 2-D input immaginary image array
%
%	Optional argument:
%	params	: an array containg 5 elements [a b c d e]
%				a = PHASE_CONSTANT
%				b = smoothing (1 for Y, 0 for N)
%				c = width of smoothing window
%				d = height of smoothing window
%				e = bandwidth of smoothing filter (variance of gaussian)
%	
%return:
%	phase	: unwrapped phase function
%
%	NTC June 26, 2008
function [A, U, V] = DiscreteDemodulation(realImg, imagImg, varargin)

optargin = size(varargin,2);
if(optargin == 1) 
	params = varargin{1};
else
	params = [300 0 5 5 0.5];
end

PHASE_CONST = params(1);
smooth = params(2);
window = [params(3) params(4)];
sigma = params(5);


% Wrap phase
wrapPhase = atan2(imagImg,realImg);


% Get image dimension
[row col] = size(realImg);

% A and wrapPhase
A = sqrt((realImg.*realImg + imagImg.*imagImg));


% Find gradient images
[realDev_H, realDev_V] = SplineGradient(realImg);
[imagDev_H, imagDev_V] = SplineGradient(imagImg);

% phase derivative (no unwrap phase needed).
Asquare = A .* A + 10^(-3);
phaseDev_H = (imagDev_H.*realImg - realDev_H.*imagImg) ./ Asquare;
phaseDev_V = (imagDev_V.*realImg - realDev_V.*imagImg) ./ Asquare;

% Smoothing
if( smooth == 1 )
	h = fspecial('gaussian', window, sigma);
	phaseDev_H = imfilter(phaseDev_H, h, 'replicate');
	phaseDev_V = imfilter(phaseDev_V, h, 'replicate');
end

U = phaseDev_H;
V = phaseDev_V;
