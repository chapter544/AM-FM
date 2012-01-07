clear all
close all
clc

%f = file2image('float', N, 'Images/barbara512R.float');
%f = file2image('float', 256, '/home/chuong/Images/lenaR.float');
%f = file2image('float', 256, '/home/chuong/Images/snake1R.float');
%f = file2image('float', 256, '/home/chuong/Images/woodR.float');
%f = file2image('float', 256, '/home/chuong/Images/burlapR.float');

%fileName = '/home/chuong/Images/EinSlack.png';
fileName = '/home/chuong/Images/lena.png';
%fileName = '/home/chuong/Images/barbara.png';
%fileName = '/home/chuong/Images/kodim05.png';
%fileName = '/home/chuong/Images/house.png';
%fileName = '/home/chuong/Images/flinstones.png';
%fileName = '/home/chuong/Images/fingerprint.png';
%fileName = '/home/chuong/Images/fingerprint.png';
f = imread(fileName, 'png');

if(size(f,3) > 1)
	f = single(f(:,:,1));
else
	f = single(f);
end

sigma = 1.0;
[u, v] = BuadesNonlinearTVCartoon(f, sigma);

figure;
subplot(2,1,1)
imagesc(u);
title(sprintf('Cartoon max=%f, min=%f\n', max(u(:)), min(u(:))));
axis image; axis off; colormap(gray(256));
subplot(2,1,2)
imagesc(v);
axis image; axis off; colormap(gray(256));
title(sprintf('Texture max=%f, min=%f\n', max(v(:)), min(v(:))));
