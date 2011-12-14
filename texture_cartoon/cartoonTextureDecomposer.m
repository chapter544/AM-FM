function [cartoon, texture, gain] = cartoonTextureDecomposer(im, alpha, beta)

%clear all
%close all
%clc

%M = 256;
%N = 256;
%M = 512;
%N = 512;
%fileName = '/home/chuong/Research/AMFM/runFall2011/Images/lenaR.float';

%clear all
%close all;
%fileName = '/home/chuong/Images/barbara512R.float';
%im = file2image2('float', 512, 512, fileName);

%fileName = '/home/chuong/Images/lenaR.float';
%im = file2image2('float', 256, 256, fileName);

%fileName = '/home/chuong/Images/fingerprint.png';
%fileName = '/home/chuong/Images/lena.png';
%fileName = '/home/chuong/Images/barbara.png';
%fileName = '/home/chuong/Images/boat.png';
%fileName = '/home/chuong/Images/peppers256.png';
%im = imread(fileName, 'png');



% convert precision to float
im = single(im);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ANALYSIS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[M,N] = size(im);
numLevels = 5;
numOrien = 6;
trans_opts.withDCRing = 1;
trans_opts.phaseUnwrap = 0;
trans_opts.PHASE_CONST = 300;
trans_opts.smoothDemodulation = 1;
trans_opts.smoothWindow = 5;
trans_opts.smoothSigma = 0.1;
trans_opts.printParams = 0;
[At Ut Vt Pt Ptls Resi] = AMFM_Transform(im, numLevels, numOrien, trans_opts);

coeff = zeros(size(im));
total = zeros(size(im));
A_total = zeros(size(im));
for levidx=1:numLevels
	AR_lev = zeros(size(im));
	for orienidx=1:numOrien
		RR = sqrt(Ut{levidx,orienidx}.^2 + Vt{levidx,orienidx}.^2);
		AR_lev = AR_lev + At{levidx,orienidx} .* RR;
		A_total = A_total + At{levidx,orienidx};
	end
	coeff(:,:,levidx) = AR_lev;
end

for levidx=1:numLevels
	coeff(:,:,levidx) = coeff(:,:,levidx) ./ A_total;
end

% calculate the gain matrix, using predefined alpha and beta
outVar = var(coeff, 1, 3);
gain = 1.0 - exp( - max(outVar-alpha, 0.0) / beta );

out = zeros(size(im));
for levidx=1:numLevels-1,
	for orienidx=1:numOrien
		out = out + At{levidx,orienidx} .* cos(Pt{levidx,orienidx});
	end
end
texture = gain .* out;
cartoon = im-texture;
%showimage(texture, 'texture');
%showimage(cartoon, 'cartoon');
