% Feed image into the cartoonTextureDecomposer and write images out to files
% function cartoonTextureDriver(rootName,fileType,M,N,inImageDir,outImageDir)
%
% Warning: no error checking on input parameter since this is just the driver
%
% Chuong nguyen
function cartoonTextureDriver(rootName, fileType, M, N, inImageDir, outImageDir)

if(nargin < 6)
	error('Not enough parameter, please type: help cartoonTextureDriver');
end

% start start timer
tstart = tic;

if(strcmp(fileType, 'float'))
	fileName = sprintf('%s/%s.float', inImageDir, rootName);
	im = file2image2('float', M, N, fileName);
	im = FullScaleStretch(im);
elseif(strcmp(fileType, 'tif'))
	fileName = sprintf('%s/%s.tif', inImageDir, rootName);
	im = imread(fileName, fileType);
elseif(strcmp(fileType, 'png'))
	fileName = sprintf('%s/%s.png', inImageDir, rootName);
	im = imread(fileName, fileType);
elseif(strcmp(fileType, 'bmp'))
	fileName = sprintf('%s/%s.bmp', inImageDir, rootName);
	im = imread(fileName, fileType);
end

% create single precision if the image has type other than 'float'
im = single(im);

dim = size(im, 3);
if(dim == 1)
	[r_cartoon, r_texture, r_gain] = cartoonTextureDecomposer(im);
	r_texture_fs = FullScaleStretch(r_texture);
	r_cartoon_fs = FullScaleStretch(r_cartoon);
	im_texture = cat(3, uint8(r_texture_fs), ...
		uint8(r_texture_fs), uint8(r_texture_fs));
	im_cartoon = cat(3, uint8(r_cartoon_fs), ...
		uint8(r_cartoon_fs), uint8(r_cartoon_fs));

	% write float type data to file
	image2file(r_texture, 'float', sprintf('%s/%s_texture.float', ...
		outImageDir, rootName), 0);
	image2file(r_cartoon, 'float', sprintf('%s/%s_cartoon.float', ...
		outImageDir, rootName), 0);
	image2file(r_gain, 'float', sprintf('%s/%s_gain.float', ...
		outImageDir, rootName), 0);
else
	[r_cartoon, r_texture, r_gain] = cartoonTextureDecomposer(im(:,:,1));
	[g_cartoon, g_texture, g_gain] = cartoonTextureDecomposer(im(:,:,2));
	[b_cartoon, b_texture, b_gain] = cartoonTextureDecomposer(im(:,:,3));

	r_texture_fs = FullScaleStretch(r_texture);
	g_texture_fs = FullScaleStretch(g_texture);
	b_texture_fs = FullScaleStretch(b_texture);

	r_cartoon_fs = FullScaleStretch(r_cartoon);
	g_cartoon_fs = FullScaleStretch(g_cartoon);
	b_cartoon_fs = FullScaleStretch(b_cartoon);

	im_texture = cat(3, uint8(r_texture_fs), ...
		uint8(g_texture_fs), uint8(b_texture_fs));
	im_cartoon = cat(3, uint8(r_cartoon_fs), ...
		uint8(g_cartoon_fs), uint8(b_cartoon_fs));

	% write to file
	image2file(r_cartoon, 'float', sprintf('%s/%s_cartoon_R.float', ...
		outImageDir, rootName), 0);
	image2file(r_texture, 'float', sprintf('%s/%s_texture_R.float', ...
		outImageDir, rootName), 0);
	image2file(r_gain, 'float', sprintf('%s/%s_gain_R.float', ...
		outImageDir, rootName), 0);

	image2file(g_cartoon, 'float', sprintf('%s/%s_cartoon_G.float', ...
		outImageDir, rootName), 0);
	image2file(g_texture, 'float', sprintf('%s/%s_texture_G.float', ...
		outImageDir, rootName), 0);
	image2file(g_gain, 'float', sprintf('%s/%s_gain_G.float', ...
		outImageDir, rootName), 0);

	image2file(b_cartoon, 'float', sprintf('%s/%s_cartoon_B.float', ...
		outImageDir, rootName), 0);
	image2file(b_texture, 'float', sprintf('%s/%s_texture_B.float', ...
		outImageDir, rootName), 0);
	image2file(b_gain, 'float', sprintf('%s/%s_gain_B.float', ...
		outImageDir, rootName), 0);
end

imwrite(im_texture, sprintf('%s/%s_texture.png', outImageDir, rootName), 'png');
imwrite(im_cartoon, sprintf('%s/%s_cartoon.png', outImageDir, rootName), 'png');

fprintf('Execution time of %dx%dx%d is %f seconds\n\n',  ...
	size(im,1), size(im,2), dim, toc(tstart));

return;
