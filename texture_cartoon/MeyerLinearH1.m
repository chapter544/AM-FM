clear all
close all
clc

%f = file2image('float', N, 'Images/barbara512R.float');
%fileName = '/home/chuong/Images/EinSlack.png';
%fileName = '/home/chuong/Images/kodim01.png';
fileName = '/home/chuong/Images/barbara.png';
f = imread(fileName, 'png');

if(size(f,3) > 1)
	f = single(f(:,:,1));
else
	f = single(f);
end

[M,N] = size(f);

% Low-pass filter L
m=(-M/2:1:M/2-1)./M;
n=(-N/2:1:N/2-1)./N;
[x,y] = meshgrid(n,m);
sigma = 1;
sigma_t = 1 / (2*pi*sigma);
L = sigma_t^4 ./ ((x.^2 + y.^2).^2 + sigma_t^4);
%L = exp( - (2*pi*sigma)^4 .* (x.^2 + y.^2).^2 );


% Compute the total variation of f_L
imDft = fftshift(fft2(f));
u = real(ifft2(ifftshift(imDft .* L)));
v = f - u;


figure;
subplot(2,2,1)
imagesc(u);
title(sprintf('Cartoon max=%f, min=%f\n', max(u(:)), min(u(:))));
axis image; axis off; colormap(gray(256));

subplot(2,2,2)
imagesc(v);
axis image; axis off; colormap(gray(256));
title(sprintf('Texture max=%f, min=%f\n', max(v(:)), min(v(:))));
subplot(2,2,3)

imagesc(f);
axis image; axis off; colormap(gray(256));
title(sprintf('Org max=%f, min=%f\n', max(f(:)), min(f(:))));

subplot(2,2,4)
imagesc(abs(L));
axis image; axis off; colormap(gray(256));
title(sprintf('L max=%f, min=%f\n', max(L(:)), min(L(:))));
