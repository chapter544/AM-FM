%phaseUnwrap2 2-D Phase Unrapping using Spline 
%Unwrap phase with derivative assistance
%paramters:
%	realImg	: 2-D input real image array
%	imagImg	: 2-D input immaginary image array
%
%	Optional argument:
%	params	: 
%				PHASE_CONST: phase constant( default = 1)
%				threshold: amplitude modulation threshold (default = 0.01) 
%	
%return:
%	phase	: unwrapped phase function
%
% 	example:
%		[A, U, V, P, Pls] = phaseUnwrap2(realImg, imagImg, PHASE_CONST, threshold)
%	NTC June 26, 2008
%	NTC August 22, 2011
%
function [A, U, V, P, Pls] = phaseUnwrap2(realImg, imagImg, varargin)

% default variables
PHASE_MULTIPLIER = 1;
A_THRESHOLD = 0.01;

optargin = size(varargin,2);
if(optargin == 1) 
	PHASE_MULTIPLIER = varargin{1};
elseif(optargin == 2)
	PHASE_MULTIPLIER = varargin{1};
	A_THRESHOLD = varargin{2};
end


% Wrap phase
wrapPhase = atan2(imagImg,realImg);


% A and wrapPhase
A = sqrt((realImg.*realImg + imagImg.*imagImg));

% Find gradient images
%[realDev_H, realDev_V] = SplineGradient(realImg);
%[imagDev_H, imagDev_V] = SplineGradient(imagImg);

[realDev_H, realDev_V] = SmoothCubicSplineGradientApproximation(realImg);
[imagDev_H, imagDev_V] = SmoothCubicSplineGradientAppriximation(imagImg);

% phase derivative (no unwrap phase needed).
Asquare = A .* A + A_THRESHOLD;
phaseDev_H = (imagDev_H.*realImg - realDev_H.*imagImg) ./ Asquare;
phaseDev_V = (imagDev_V.*realImg - realDev_V.*imagImg) ./ Asquare;


% least-square phase unwrapping algorithm
unWrapPhase = poisson_solver_function_neumann(phaseDev_H, phaseDev_V);
Pls = unWrapPhase;

% Find congruence unwrapped phase
B = zeros(size(realImg));
B = floor((PHASE_MULTIPLIER.* unWrapPhase - wrapPhase) ./ (2*pi));
P = wrapPhase + (2*pi) .* floor(B);

% Compute the gradient
[U,V] = FirstOrderGradient(P);

U = U ./ PHASE_MULTIPLIER;
V = V ./ PHASE_MULTIPLIER;
