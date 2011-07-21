%phaseUnwrap2 2-D Phase Unrapping using Spline 
%Unwrap phase with derivative assistance
%paramters:
%	realImg	: 2-D input real image array
%	imagImg	: 2-D input immaginary image array
%
%	Optional argument:
%	params	: an array containg 6 elements [a b c d e f]
%				a = PHASE_CONSTANT
%				b = smoothing (1 for Y, 0 for N)
%				c = width of smoothing window
%				d = height of smoothing window
%				e = bandwidth of R (radial frequency)
%				f = bandwidth of T (angular frequency)
%	
%return:
%	phase	: unwrapped phase function
%
%	NTC June 26, 2008
function [A, U, V, phase, unWrapPhase] = phaseUnwrapAdjusted(realImg, imagImg, varargin)

optargin = size(varargin,2);
if(optargin == 1) 
	params = varargin{1};
else
	params = [1 0 5 5 1 1];
end

PHASE_CONST = params(1);
smooth = params(2);
window = [params(3) params(4)];
alpha = params(5);
beta = params(6);


% Wrap phase
wrapPhase = atan2(imagImg,realImg);


LOW_THRESHOLD = 10^(-5);

% Get image dimension
[row col] = size(realImg);

% A and wrapPhase
A = sqrt((realImg.*realImg + imagImg.*imagImg));
Asquare = A .* A + LOW_THRESHOLD .* (A < LOW_THRESHOLD); 


% Find gradient images
[realDev_H, realDev_V] = SplineGradient(realImg);
[imagDev_H, imagDev_V] = SplineGradient(imagImg);


% phase derivative (no unwrap phase needed).
phaseDev_H = (imagDev_H.*realImg - realDev_H.*imagImg) ./ Asquare;
phaseDev_V = (imagDev_V.*realImg - realDev_V.*imagImg) ./ Asquare;


% least-square phase unwrapping algorithm
unWrapPhase = poisson_solver_function_neumann(phaseDev_H, phaseDev_V);
showimage(unWrapPhase);

%showimage(unWrapPhase)

% Find congruence unwrapped phase
B = zeros(row,col);
B = floor((PHASE_CONST .* unWrapPhase - wrapPhase) ./ (2*pi));
phase = wrapPhase + (2*pi) .* floor(B);


%H = fspecial('gaussian', window, 1.0);
%Bfilt = imfilter(B, H, 'replicate');
%phase = wrapPhase + (2*pi) .* floor(Bfilt);

% Compute the gradient
[U,V] = FirstOrderGradient(phase);

U = U ./ PHASE_CONST;
V = V ./ PHASE_CONST;
