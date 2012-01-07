function [u, v] = BuadesNonlinearTVCartoon(f, sigma)

[M,N] = size(f);

% Low-pass filter L
m=(-M/2:1:M/2-1)./M;
n=(-N/2:1:N/2-1)./N;
[x,y] = meshgrid(n,m);
%sigma = 1;
%sigma_t = 1 / (2*pi*sigma);
r_square = x.^2 + y.^2;
%L = sigma_t^4 ./ (r_square.^2 + sigma_t^4);
L = exp( - (2*pi*sigma)^4 .* r_square.^2 );

% Gaussian filter with std sigma
g = fspecial('gaussian', [7 7], sigma);

% Compute total variation of f
[f_x, f_y] = gradient(f);
f_TV = sqrt(f_x.^2 + f_y.^2);
f_TV_G = imfilter(f_TV, g, 'replicate');

% Compute the total variation of f_L
imDft = fftshift(fft2(f));
f_L = real(ifft2(ifftshift(imDft .* L)));
[f_L_x, f_L_y] = gradient(f_L);
f_L_TV = sqrt(f_L_x.^2 + f_L_y.^2);
f_L_TV_G = imfilter(f_L_TV, g, 'replicate');

% Compute lamda
lamda = (f_TV_G - f_L_TV_G) ./ f_TV_G;

showimage(lamda)

u = tv_weight(lamda) .* f_L + (1 - tv_weight(lamda)) .* f;
v = f - u;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function y = tv_weight(x)
[M,N] = size(x);
y = zeros(size(x));

a1 = 0.25;
a2 = 0.5;

for m=1:M,
	for n=1:N,
		if(x(m,n) < a1);
				y(m,n) = 0.0;
		elseif( x(m,n) > a2);
				y(m,n) = 1.0;
		else
				y(m,n) = 1 / (a2 - a1) * x(m,n) - 1;
		end
	end
end
end


