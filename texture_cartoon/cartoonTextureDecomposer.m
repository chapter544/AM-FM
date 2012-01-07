% Decompose an image into texture + cartoon
% function [cartoon, texture, gain] = cartoonTextureDecompoer
%
% input:
%	im:	2D input image
%
% outputs:
%	cartoon: 2D cartoon image
%	texture: 2D texture image
%	gain: gain mask with range [0,1]
%
% NTC December 2011
function [cartoon, texture, gain] = cartoonTextureDecomposer(im)

% default paramters
% wL = 7 : smoothing window for gradient computation
% THRESHOLD: 0.25 --> hardcoded THRESHOLD
% standard trans_opts: with 5 levels and 6 orientations

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


% compute the estimation of image gradient
dx = zeros(size(im));
dy = zeros(size(im));
dx_t = zeros(size(im));
dy_t = zeros(size(im));
for levidx=1:numLevels
	dx_lev = zeros(size(im));
	dy_lev = zeros(size(im));
	for orienidx=1:numOrien
		dx_lev = dx_lev + Ut{levidx,orienidx} .* At{levidx,orienidx} .* sin(Pt{levidx,orienidx});
		dy_lev = dy_lev + Vt{levidx,orienidx} .* At{levidx,orienidx} .* sin(Pt{levidx,orienidx});
	end

	dx(:,:,levidx) = dx_lev;
	dy(:,:,levidx) = dy_lev;
	dx_t = dx_t + dx_lev;
	dy_t = dy_t + dy_lev;
end


% smoothing out the total gradient
wL = 7;
h = fspecial('gaussian', [wL wL], 3.0);
totalGrad0 = sqrt(dx_t.^2 + dy_t.^2);
totalGrad = imfilter(totalGrad0, h, 'replicate');

% gradient estimation per level
levGrad = zeros(size(im));
for idx=1:numLevels
	xx = zeros(size(im));
	yy = zeros(size(im));
	for levidx=idx:numLevels
		xx = xx + dx(:,:,levidx);	
		yy = yy + dy(:,:,levidx);	
	end
	r = sqrt(xx.^2 + yy.^2);
	r_filt = imfilter(r, h, 'replicate');
	levGrad(:,:,idx) = (totalGrad - r_filt) ./ totalGrad;
end


% get the masking level, the threshold is set to 0.25 hard
THRESHOLD = 0.25;
mask = zeros(size(im));
for m=1:M
	for n=1:N
		t = levGrad(m,n,1);
		for levidx=1:numLevels
			t = t + levGrad(m,n,levidx);
			if(t > THRESHOLD)
				mask(m,n) = levidx;
				break;
			end
		end
	end
end

% the gain is computed by first finding the mask
myChan = median(mask(:));
gain = 0.0*(levGrad(:,:,myChan) < THRESHOLD) + 1.0 * (levGrad(:,:,myChan) >= THRESHOLD);
gain = imfilter(gain, h, 'replicate');

% compute the final textural output
out = zeros(size(im));
for m=1:M
	for n=1:N
		myChan = mask(m,n);
		if(myChan >= numLevels)
			myChan = -1;
		end
		for levidx=1:myChan
			for orienidx=1:numOrien
				out(m,n) = out(m,n) + At{levidx,orienidx}(m,n) .* cos(Pt{levidx,orienidx}(m,n));
			end
		end
	end
end
texture = out .* gain;
cartoon = im-texture;

%subplot(2,2,1)
%imagesc(im); axis('image'); axis off; colormap(gray(256));
%subplot(2,2,2)
%imagesc(gain); axis('image'); axis off; colormap(gray(256));
%subplot(2,2,3)
%imagesc(texture); axis('image'); axis off; colormap(gray(256));title('texture');
%subplot(2,2,4)
%imagesc(cartoon); axis('image'); axis off; colormap(gray(256));title('cartoon');
