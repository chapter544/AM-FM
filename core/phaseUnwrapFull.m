function [A,U,V,P,Pls] = phaseDemodFull(realImg, imagImg, varargin)

optargin = size(varargin,2);
if(optargin == 1) 
	params = varargin{1};
else
	params = [0 5 5 1];
end

smooth = params(1);
window = [params(2) params(3)];
sigma = params(4);

% Wrap phase
wrapPhase = atan2(imagImg,realImg);

LOW_THRESHOLD = 10^(-6);

% Get image dimension
[M N] = size(realImg);

% A and wrapPhase
A = sqrt((realImg.*realImg + imagImg.*imagImg));

% Find gradient images
[realDev_H, realDev_V] = SplineGradient(realImg);
[imagDev_H, imagDev_V] = SplineGradient(imagImg);

% phase derivative (no unwrap phase needed).
Asquare = A .* A + LOW_THRESHOLD;
phaseDev_H = (imagDev_H.*realImg - realDev_H.*imagImg) ./ Asquare;
phaseDev_V = (imagDev_V.*realImg - realDev_V.*imagImg) ./ Asquare;


% Smoothing
if( smooth == 1 )
	h = fspecial('gaussian', window, sigma);
	phaseDev_H = imfilter(phaseDev_H, h, 'replicate');
	phaseDev_V = imfilter(phaseDev_V, h, 'replicate');
end

% least-square phase unwrapping algorithm
unWrapPhase = poisson_solver_function_neumann(phaseDev_H, phaseDev_V);

% Find congruence unwrapped phase
B = zeros(M,N);
B = floor((unWrapPhase - wrapPhase) ./ (2*pi));
P = wrapPhase + (2*pi) .* floor(B);
[U,V] = FirstOrderGradient(P);

Pls = unWrapPhase;
