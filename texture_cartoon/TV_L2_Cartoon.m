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
m=(-M/2:1:M/2-1)/M;
n=(-N/2:1:N/2-1)/N;
[x,y] = meshgrid(n,m);
sigma = 4;
sigma_t = 1 / (2*pi*sigma);
L = sigma_t^4 ./ ((x.^2 + y.^2).^2 + sigma_t^4);

b = 1;
a = 0;
m = min(f(:));
M = max(f(:));
f = (b-a) * (f-m)/(M-m) + a;

options.dt = 0.1249;
lam = 1/20 * 215 ./ max(f(:));
options.niter = 300;
%options.use_gabor = -1; 
% 1:TV-Hilbert, -1:TV-L2
options.use_gabor = -1; 
options.p0x = []; options.p0y = [];
disp('--> Performing TV-L2 decomposition.');
options.p0x = []; options.p0y = [];
[u0,v0,options.p0x,options.p0y] = perform_tv_hilbert_projection(f,[],lam,options);

image2file(u0, 'float', sprintf('%s/%s_cartoon.png', ...
	outImageDir, rootName), 0);
image2file(r_cartoon, 'float', sprintf('%s/%s_texture.png', ...
	outImageDir, rootName), 0);
image2file(r_gain, 'float', sprintf('%s/%s_gain.png', ...
	outImageDir, rootName), 0);




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

%u = zeros(size(f));
%lambda = 0.5;
%K = 30;
%for idx=1:K
%	imDft = fftshift(fft2(f - u));
%	[ux, uy] = gradient(u);
%	mag_grad = sqrt(ux.^2 + uy.^2);
%	mag_grad = mag_grad + (mag_grad == 0);
%	ux = ux ./ mag_grad;
%	uy = uy ./ mag_grad;
%	[uxx, uxy] = gradient(ux);
%	[uyx, uyy] = gradient(uy);
%
%	firstterm = uxx + uyy;
%	%secondterm = real(ifft2(ifftshift(imDft .* L .* L)));
%	secondterm = real(ifft2(ifftshift(imDft .* L)));
%	%showimage(firstterm);
%	%showimage(secondterm);
%	u = firstterm +  lambda * secondterm;
%end
%
%v = f - u;
%
%
%
%b = 1;
%a = 0;
%m = min(f(:));
%M = max(f(:));
%f = (b-a) * (f-m)/(M-m) + a;
%
%
