function [A Uest Vest P Pls] = RieszPhaseUnwrap(f, varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [A U V P Pls]  = RieszPhaseUnwrap(f, PHASE_CONST, smoothing, w)
%	f			: 2-D input function
%	PHASE_CONST	: phase constant (default = 300)
%	smoothing	: use smooth filter (value = {0, 1}, default = 1)
%	w			: width of the smoothing kernel
%	
%	A		: amplitude modulation function (AM)
%	Uest	: horizontal frequency function (horizontal FM component)
%	Vest	: vertical frequency function (vertical FM component)
%	P		: true phase function
%	Pls		: Least square phase function
%
%	NTC 09/27/2010
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if(nargin < 1)
	error('Missing input image');
end


% pre-defined parameters
epsilon = 10^(-3); % add constant to denominator
PHASE_CONST = 300.0;
smoothing = 1;
w = 3;

optargin = size(varargin,2);
if(optargin == 1) 
	PHASE_CONST = varargin{1};
end
if( optargin == 2 )
	smoothing = varargin{2};
	if(smoothing ~= 1)
		smoothing = 0;
	end
end

if( optargin == 3 )
	w = varargin{3};
	if(w < 3)
		w = 3;
	end
end


% size of input
%[M, N] = size(f);

% padded f
L = 16;
f = padarray(f, [L L], 'symmetric', 'both');
[M, N] = size(f);

% create two riez filters
Hu = zeros(M,N);
Hv = zeros(M,N);
for m=-M/2:1:M/2-1
	for n=-N/2:1:N/2-1
		if(n == 0 && m == 0)
			Hu(M/2+1,N/2+1) = 0.0;
			Hv(M/2+1,N/2+1) = 0.0;
		else
			Hu(m+M/2+1,n+N/2+1) = n / sqrt(m^2 + n^2);
			Hv(m+M/2+1,n+N/2+1) = m / sqrt(m^2 + n^2);
		end
	end
end

% a is the horizontal component of the Riesz transform
f_a = imag(ifft2(ifftshift(fftshift(fft2(f)) .* Hu)));
% b is the vertical component of the Riesz transform
f_b = imag(ifft2(ifftshift(fftshift(fft2(f)) .* Hv)));

% Riesz AM function, isotropic kernel
A = sqrt(f.^2 + f_a.^2 + f_b.^2);

% absolute value of the imaginary image
f_im_abs = sqrt(f_a.^2 + f_b.^2);


% Calculate the magnitude and argument of the FM function
% using the demodulation algorithm applying straight tounWrapPhasen 
% the contitunous equation
%[Ax, Ay] = SplineGradient(A);
%[fx, fy] = SplineGradient(f);

% The following use Farid and Simoncelli derivative filter
[fx, fy] = derivative5(f, 'x', 'y');
[Ax, Ay] = derivative5(A, 'x', 'y');

deno = A .* f_im_abs;
%if(smoothing == 1)
%	h = fspecial('gaussian', [w w], 0.5);
%	deno = imfilter(deno, h, 'replicate');
%end
Uabs = abs(Ax .* f - A .* fx) ./ (deno + epsilon);
Vabs = abs(Ay .* f - A .* fy) ./ (deno + epsilon);

%Uabs = abs(Ax .* f - A .* fx) ./ (A .* f_im_abs + epsilon);
%Vabs = abs(Ay .* f - A .* fy) ./ (A .* f_im_abs + epsilon);

% get rid of some outliers by smoothing
if( smoothing == 1)
	Uabs = medfilt2(Uabs, [w w]);
	Vabs = medfilt2(Vabs, [w w]);
end
R = sqrt(Uabs.^2 + Vabs.^2);

% find the angurment of the FM, limit to (-pi/2, pi/2].
theta = atan( (Ay .* f - A .* fy) ./  (Ax .* f - A .* fx + epsilon) );
for m=1:M
	for n=1:N
		if( abs(theta(m,n) + pi/2) < 0.01 )
			theta(m,n) = theta(m,n) + pi;
		end
	end
end
if( smoothing == 1)
	theta = medfilt2(theta, [w w]);
end



% Perform least square phase intergration 
Uest = R .* cos(theta);
Vest = R .* sin(theta);
Pls = poisson_solver_function_neumann(Uest, Vest);


% enforce phase congruency for true phase function
% wrapPhase = real(acos(f ./ (A + epsilon)));
% NOTE: change from using acos to using atan2 with f_im_abs
wrapPhase = atan2(f_im_abs, f);
B = zeros(M,N);
B = floor((PHASE_CONST .* Pls - wrapPhase) ./ (2*pi));
P = wrapPhase + (2*pi) .* floor(B);

% get the un-padded output
Uest = Uest(L+1:M-L, L+1:N-L);
Vest = Vest(L+1:M-L, L+1:N-L);
A = A(L+1:M-L, L+1:N-L);
P = P(L+1:M-L, L+1:N-L);
Pls = Pls(L+1:M-L, L+1:N-L);
