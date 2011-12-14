function cartoonTextureDriver(rootName, fileType, M, N, alpha, beta, inImageDir, outImageDir)

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

im = single(im);
dim = size(im, 3);
if(dim == 1)
	[r_cartoon, r_texture, r_gain] = cartoonTextureDecomposer(im, alpha, beta);
	r_texture_fs = FullScaleStretch(r_texture);
	r_cartoon_fs = FullScaleStretch(r_cartoon);
	im_texture = cat(3, uint8(r_texture_fs), ...
		uint8(r_texture_fs), uint8(r_texture_fs));
	im_cartoon = cat(3, uint8(r_cartoon_fs), ...
		uint8(r_cartoon_fs), uint8(r_cartoon_fs));

	% write float type data to file
	image2file(r_cartoon, 'float', sprintf('%s/%s_texture.png', ...
		outImageDir, rootName), 0);
	image2file(r_cartoon, 'float', sprintf('%s/%s_texture.png', ...
		outImageDir, rootName), 0);
	image2file(r_gain, 'float', sprintf('%s/%s_gain.png', ...
		outImageDir, rootName), 0);
else
	[r_cartoon, r_texture, r_gain] = cartoonTextureDecomposer(im(:,:,1), alpha, beta);
	[g_cartoon, g_texture, g_gain] = cartoonTextureDecomposer(im(:,:,2), alpha, beta);
	[b_cartoon, b_texture, b_gain] = cartoonTextureDecomposer(im(:,:,3), alpha, beta);

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
	image2file(r_cartoon, 'float', sprintf('%s/%s_texture_R.png', ...
		outImageDir, rootName), 0);
	image2file(r_cartoon, 'float', sprintf('%s/%s_texture_R.png', ...
		outImageDir, rootName), 0);
	image2file(r_gain, 'float', sprintf('%s/%s_gain_R.png', ...
		outImageDir, rootName), 0);

	image2file(g_cartoon, 'float', sprintf('%s/%s_texture_G.png', ...
		outImageDir, rootName), 0);
	image2file(g_cartoon, 'float', sprintf('%s/%s_texture_G.png', ...
		outImageDir, rootName), 0);
	image2file(g_gain, 'float', sprintf('%s/%s_gain_G.png', ...
		outImageDir, rootName), 0);

	image2file(b_cartoon, 'float', sprintf('%s/%s_texture_B.png', ...
		outImageDir, rootName), 0);
	image2file(b_cartoon, 'float', sprintf('%s/%s_texture_B.png', ...
		outImageDir, rootName), 0);
	image2file(b_gain, 'float', sprintf('%s/%s_gain_B.png', ...
		outImageDir, rootName), 0);
end

imwrite(im_texture, sprintf('%s/%s_texture.png', outImageDir, rootName), 'png');
imwrite(im_cartoon, sprintf('%s/%s_cartoon.png', outImageDir, rootName), 'png');

fprintf('Execution time of %dx%dx%d is %f seconds\n\n',  M, N, dim, toc(tstart));

return;
