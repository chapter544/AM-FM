%phaseUnwrap2 2-D Phase Unrapping using Spline 
%Unwrap phase with derivative assistance
%paramters:
%	realImg	: 2-D input real image array
%	imagImg	: 2-D input immaginary image array
%
%	Optional argument:
%	params	: an array containg 5 elements [a b c d e]
%				a = perform phase unwrapping
%				b = PHASE_CONSTANT
%				c = smoothing (1 for Y, 0 for N)
%				d = size of smoothing window
%				e = bandwidth of smoothing filter (variance of gaussian)
%	
%return:
%	phase	: unwrapped phase function
%
%	NTC June 26, 2008
%	NTC December 2, 2011
function [A, U, V, P, Pls] = phaseUnwrap(realImg, imagImg, varargin)
if( size(realImg,1) ~= size(imagImg,1) || ...
    size(realImg,2) ~= size(imagImg,2) )
    error('Dimension of realImg and imagImg mismatch.');
end

% demodulation parameters
optargin = size(varargin,2);
if(optargin == 1) 
	params = varargin{1};
else
	params = [1 300 1 3 0.5];
end

performPhaseUnwrapping = params(1);
PHASE_CONST = params(2);
smoothing = params(3);
window = [params(4) params(4)];
sigma = params(5);


% outputs allocations
A = zeros(size(realImg));
U = zeros(size(realImg));
V = zeros(size(realImg));
P = zeros(size(realImg));
Pls = zeros(size(realImg));


%fprintf('Performing AM-FM demodulation using:\n');
%fprintf('\tPHASE_CONST: %0.2f\n',PHASE_CONST);
%fprintf('\tsmoothing: %d\n',smooth);
%if(smooth == 1) 
%	fprintf('\nKernel: Gaussian\n');
%	fprintf('\t\twindow: %dx%d\n',window(1), window(2));
%	fprintf('\t\tsigma: %0.2f\n', sigma);
%end

% A and wrapPhase
A = sqrt((realImg.*realImg + imagImg.*imagImg));

% Find gradient images
[realDev_H, realDev_V] = SplineGradient(realImg);
[imagDev_H, imagDev_V] = SplineGradient(imagImg);

% phase derivative (no unwrap phase needed).
THRESHOLD = 10^(-3);
Asquare = A .* A;
Asquare = Asquare + (Asquare < THRESHOLD) .* THRESHOLD;
phaseDev_H = (imagDev_H.*realImg - realDev_H.*imagImg) ./ Asquare;
phaseDev_V = (imagDev_V.*realImg - realDev_V.*imagImg) ./ Asquare;

% Smoothing
if( smoothing == 1 )
	h = fspecial('gaussian', window, sigma);
	phaseDev_H = imfilter(phaseDev_H, h, 'replicate');
	phaseDev_V = imfilter(phaseDev_V, h, 'replicate');
elseif( smoothing == 2)
    phaseDev_H = medfilt2(phaseDev_H, window);
	phaseDev_V = medfilt2(phaseDev_V, window);
end


% compute  wrapped phase from the real and imaginary parts
wrapPhase = atan2(imagImg,realImg);
if(performPhaseUnwrapping == 1)
	% least-square phase unwrapping algorithm
	Pls = poisson_solver_function_neumann(phaseDev_H, phaseDev_V);

	% Find congruence unwrapped phase, by computing PHASE_CONST * Pactual
	B = floor((PHASE_CONST .* Pls - wrapPhase) ./ (2*pi));
	P = wrapPhase + (2*pi) .* floor(B);

	% Compute the gradient with scaled-back P
	[U,V] = FirstOrderGradient(P./PHASE_CONST);

	% scale phase back
	%U = U ./ PHASE_CONST;
	%V = V ./ PHASE_CONST;
else
	U = phaseDev_H;
	V = phaseDev_V;
	P = wrapPhase;
end
