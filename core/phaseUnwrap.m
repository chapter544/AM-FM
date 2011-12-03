% AM-FM demodulation 
%	Compute the AM and FM functions from the real and imaginary image.
%
% function [A, U, V, P, Pls] = phaseUnwrap(realImg, imagImg, varargin)
% paramters:
%	realImg	: 2-D input real image array
%	imagImg	: 2-D input immaginary image array
%
%	Optional argument: demod_opts
%	demod_opts.phaseUnwrap = 1;  --> perform phase unwrapping
%	demod_opts.PHASE_CONST = 300;  --> phase scaling for smooth phase
%	demod_opts.smoothDemodulation = 1; --> smoothing estimated FM functions
%	demod_opts.smoothWindow = 3; --> smoothing window 
%	demod_opts.smoothSigma = 0.5; --> smoothing sigma (Gaussian kernel)
% return:
%	phase	: unwrapped phase function
%
%	NTC June 26, 2008
%	NTC December 2, 2011
%
function [A, U, V, P, Pls] = phaseUnwrap(realImg, imagImg, varargin)
if( size(realImg,1) ~= size(imagImg,1) || ...
    size(realImg,2) ~= size(imagImg,2) )
    error('Dimension of realImg and imagImg mismatch.');
end

% demodulation parameters
optargin = size(varargin,2);
if(optargin == 1) 
	demod_opts = varargin{1};
else
	demod_opts.phaseUnwrap = 1; 
	demod_opts.PHASE_CONST = 300; 
	demod_opts.smoothDemodulation = 1;
	demod_opts.smoothWindow = 3;
	demod_opts.smoothSigma = 0.5;
	demod_opts.printParams = 1;
end

% outputs allocations
A = zeros(size(realImg));
U = zeros(size(realImg));
V = zeros(size(realImg));
P = zeros(size(realImg));
Pls = zeros(size(realImg));

if(demod_opts.printParams == 1)
	fprintf('Performing AM-FM demodulation using:\n');
	fprintf('\tPHASE_CONST: %0.2f\n',demod_opts.PHASE_CONST);
	fprintf('\tsmoothing demodulation: %d\n',demod_opts.smoothDemodulation);
	if(demod_opts.smoothDemodulation == 1) 
		fprintf('\t\tKernel: Gaussian\n');
		fprintf('\t\twindow: %dx%d\n', demod_opts.smoothWindow, ...
			demod_opts.smoothWindow);
		fprintf('\t\tsigma: %0.2f\n', demod_opts.smoothSigma);
	end
end


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
smoothingWindow = [demod_opts.smoothWindow demod_opts.smoothWindow];
if( demod_opts.smoothDemodulation == 1 ) % gaussian kernel
	h = fspecial('gaussian', smoothingWindow, demod_opts.smoothSigma);
	phaseDev_H = imfilter(phaseDev_H, h, 'replicate');
	phaseDev_V = imfilter(phaseDev_V, h, 'replicate');
elseif( demod_opts.smoothingDemodulation == 2) % median filter
    phaseDev_H = medfilt2(phaseDev_H, smoothingWindow);
	phaseDev_V = medfilt2(phaseDev_V, smoothingWindow);
end


% compute  wrapped phase from the real and imaginary parts
wrapPhase = atan2(imagImg,realImg);
if(demod_opts.phaseUnwrap == 1)
	% least-square phase unwrapping algorithm
	Pls = poisson_solver_function_neumann(phaseDev_H, phaseDev_V);

	% Find congruence unwrapped phase, by computing PHASE_CONST * Pactual
	B = floor((demod_opts.PHASE_CONST .* Pls - wrapPhase) ./ (2*pi));
	P = wrapPhase + (2*pi) .* floor(B);

	% Compute the gradient with scaled-back P
	[U,V] = FirstOrderGradient(P./demod_opts.PHASE_CONST);
	% scale phase back
	%U = U ./ PHASE_CONST;
	%V = V ./ PHASE_CONST;
else
	U = phaseDev_H;
	V = phaseDev_V;
	P = wrapPhase;
end
