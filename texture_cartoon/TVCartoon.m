clear all
close all
clc

%f = file2image('float', N, 'Images/barbara512R.float');
%fileName = '/home/chuong/Images/EinSlack.png';
fileName = '/home/chuong/Images/kodim01.png';
f = imread(fileName, 'png');

if(size(f,3) > 1)
	f = single(f(:,:,1));
else
	f = single(f);
end

[M,N] = size(f);

% Low-pass filter L
m=-M/2:1:M/2-1;
n=-N/2:1:N/2-1;
[x,y] = meshgrid(n,m);
sigma = 16;
L = sigma^4 ./ ((x.^2 + y.^2).^2 + sigma^4);

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
u = tv_weight(lamda) .* f_L + (1 - tv_weight(lamda)) .* f;
v = f - u;

figure;
subplot(2,1,1)
imagesc(u);
title(sprintf('Cartoon max=%f, min=%f\n', max(u(:)), min(u(:))));
axis image; axis off; colormap(gray(256));
subplot(2,1,2)
imagesc(v);
axis image; axis off; colormap(gray(256));
title(sprintf('Texture max=%f, min=%f\n', max(v(:)), min(v(:))));

