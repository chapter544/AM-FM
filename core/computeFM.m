% Compute the A, U, V without phase unwrapping
function [A, U, V] = computeFM(realImg, imagImg, varargin)

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
U = (imagDev_H.*realImg - realDev_H.*imagImg) ./ Asquare;
V = (imagDev_V.*realImg - realDev_V.*imagImg) ./ Asquare;


