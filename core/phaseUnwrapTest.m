%phaseUnwrap2 2-D Phase Unrapping using Spline 
%Unwrap phase with derivative assistance
%paramters:
%	realImg	: 2-D input real image array
%	imagImg	: 2-D input immaginary image array
%	
%return:
%	phase	: unwrapped phase function
%
%	NTC June 26, 2008
function [A, U, V, phase] = phaseUnwrap(realImg, imagImg, params)

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

%imagImgCS_R = ImageCubicSplineOnDirection(imagImg,1);
%imagImgCS_C = ImageCubicSplineOnDirection(imagImg,2);

% Find gradient images
[realDev_H, realDev_V] = SplineGradient(realImg);
[imagDev_H, imagDev_V] = SplineGradient(imagImg);

% phase derivative (no unwrap phase needed).
phaseDev_H = (imagDev_H.*realImg - realDev_H.*imagImg) ./ Asquare;
phaseDev_V = (imagDev_V.*realImg - realDev_V.*imagImg) ./ Asquare;

showimage(realImg, 'real');
showimage(imagImg, 'imag');
showimage(phaseDev_H)
showimage(phaseDev_V)

m_needle(phaseDev_H, phaseDev_V, [], 256, -1, '', 8, [])



%R_est = sqrt(phaseDev_H.^2 + phaseDev_V.^2);
%T_est = atan2(phaseDev_V, phaseDev_H);

%% Post processing
%if(smooth == 1) 
%	std_R = std(R_est(:));
%	mean_R = mean(R_est(:));
%	std_T = std(T_est(:));
%	mean_T = mean(T_est(:));
%
%	for m=2:row
%		for n=2:col
%			if(abs(R_est(m,n) - mean_R) > alpha.*std_R) 
%				up = R_est(m-1,n);
%				left = R_est(m,n-1);
%				top_left = R_est(m-1,n-1);
%				R_est(m,n) = (up + left + top_left) / 3.0;
%			end
%
%			if(abs(T_est(m,n) - mean_T) > beta.*std_T) 
%				up = T_est(m-1,n);
%				left = T_est(m,n-1);
%				top_left = T_est(m-1,n-1);
%				T_est(m,n) = (up + left + top_left) / 3.0;
%			end
%		end % end for
%	end % end for
%end % end if
%phaseDev_H = medfilt2(phaseDev_H, window);
%phaseDev_V = medfilt2(phaseDev_V, window);

%phaseDev_H = R_est .* cos(T_est);
%phaseDev_V = R_est .* sin(T_est);
%phaseDev_H = abs(phaseDev_H);

% least-square phase unwrapping algorithm
unWrapPhase = poisson_solver_function_neumann(phaseDev_H, phaseDev_V);

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
