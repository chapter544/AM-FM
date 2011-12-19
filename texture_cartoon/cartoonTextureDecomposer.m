% Cartoon Texture Decomposition
% function [cartoon, texture, gain] = cartoonTextureDecomper(im)
% inputs:
%		im: 2-D image, limited to 1 channel for now
%
% outputs:
%		cartoon: cartoon component of image
%		texture: texture component of image
%		gain: the gain function [0,1], higher gain means texture  
%
% Chuong Nguyen
%
function [cartoon, texture, gain] = cartoonTextureDecomposer(im)


if(nargin < 1)
	error('Not enough parameter!');
end

% for now, just do 1 channel at a time (2011/12/18)
if(size(im,3) > 1)
	im = im(:,:,1);
end

% predefine a [3x3] smoothing window to fix NaN numbers
h = fspecial('average', [3 3]);

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

coeffR = zeros(size(im));
A_total = zeros(size(im));
L = numLevels-1;
for levidx=1:L
	AR_lev = zeros(size(im));
	for orienidx=1:numOrien
		RR = sqrt(Ut{levidx,orienidx}.^2 + Vt{levidx,orienidx}.^2);
		AR_lev = AR_lev + At{levidx,orienidx} .* RR;
		A_total = A_total + At{levidx,orienidx};
	end
	coeffR(:,:,levidx) = AR_lev;
end

% normalize with A_total
for levidx=1:L
	coeffR(:,:,levidx) = coeffR(:,:,levidx) ./ A_total;
end


% compute variance
outVar = var(coeffR, 1, 3);
outVar(isnan(outVar)) = 0.0;


% compute weighted mean for starting texture levels
coeffRw = zeros(size(im));
for levidx=1:L
	coeffRw(:,:,levidx) = coeffR(:,:,levidx) .* levidx;
end
outMean = sum(coeffRw, 3) ./ sum(coeffR,3);
outMean(isnan(outMean)) = 0.0;
% smooth it out so that isnan is fixed
filterChan  = imfilter(outMean, h, 'symmetric');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compute gain 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% estimate alpha using the Rayliegh distribution

alpha = sqrt(0.5*pi) * sqrt(sum(sum(outVar.^2)) / (2*prod(size(outVar))));
gain0 = exp( -max(alpha - outVar, 0.0) / (0.5*alpha) );
% weight the e
gain = imfilter(gain0, h, 'symmetric');

% extracting texture and cartoon
out = zeros(size(im));
for m=1:M
	for n=1:N
		chan = ceil(filterChan(m,n));
		for levidx=1:chan
			for orienidx=1:numOrien
				out(m,n) = out(m,n) + ...
					At{levidx,orienidx}(m,n) .* cos(Pt{levidx,orienidx}(m,n));
			end
		end
	end
end
texture = gain .* out;
cartoon = im-texture;
%showimage(texture, 'texture');
%showimage(cartoon, 'cartoon');
