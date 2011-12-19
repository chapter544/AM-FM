clear all
close all
clc

%f = file2image('float', N, 'Images/barbara512R.float');
%fileName = '/home/chuong/Images/EinSlack.png';
fileName = '/home/chuong/Images/barbara.png';
%fileName = '/home/chuong/Images/kodim01.png';
f = imread(fileName, 'png');

if(size(f,3) > 1)
	f = single(f(:,:,1));
else
	f = single(f);
end

sigma = 8;
[u, v] = NonlinearTVCartoon(f, sigma);

figure;
subplot(2,1,1)
imagesc(u);
title(sprintf('Cartoon max=%f, min=%f\n', max(u(:)), min(u(:))));
axis image; axis off; colormap(gray(256));
subplot(2,1,2)
imagesc(v);
axis image; axis off; colormap(gray(256));
title(sprintf('Texture max=%f, min=%f\n', max(v(:)), min(v(:))));
