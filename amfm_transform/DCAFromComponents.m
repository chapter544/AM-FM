% Perform Dominant component analysis from AM-FM components
% Perform pixel-wise dominant component with maximizing criteria A*R
%
% function [domA domU domV domP] = DCAFromComponents(A, U, V, P, dca_opts)
%
%	inputs:
%			A: structure cell ( numLevels x numOrien x size(im) )
%			U: structure cell ( numLevels x numOrien x size(im) )
%			V: structure cell ( numLevels x numOrien x size(im) )
%			P: structure cell ( numLevels x numOrien x size(im) )
%
%	dca_opts is an object:
%	dca_opts.useNeighbor = 1 for weighted average of neighbor pixels, 
%						   0 is NO.  
%	dca_opts.window = 3 --> size of neighbor window.  
%	dca_opts.useWeightedAMFM = 1 --> the dominant pixel
%		is the product of AM and FM, 0 for AM only.
%	dca_opts.K --> number of levels for DCA computation
%
%
%	outputs:
%			domA:	dominant AM
%			domU:	dominant FM horizontal direction
%			domV:	dominant FM horizontal direction
%			domP:	dominant phase function
%			idx:	array containing dominant orientation
%
function [domA domU domV domP idx] = DCAFromComponents(...
									strutA, strutU, strutV, strutP, varargin);

if( nargin < 4 )
	error('Not enough input argument: DCAFromComponents(A, U, V, P, dca_opts)');
end

[numLevels, numOrien] = size(strutA);

% Get parameters
if( length(varargin) == 1 )
	dca_opts = varargin{1};
else
	dca_opts.useNeighbor = 1;
	dca_opts.window = 3;
	dca_opts.useWeightedAMFM = 1;
	dca_opts.K = numLevels;
end
if(dca_opts.K <= 0 || dca_opts.K > numLevels)
	dca_opts.K = numLevels;
end


% Regenerate these structures into 3-D array
[M,N] = size(strutA{1,1});
A = zeros(M,N);
U = zeros(M,N);
V = zeros(M,N);
P = zeros(M,N);
R = zeros(M,N);
for levidx=1:dca_opts.K,
	for orienidx=1:numOrien,
		A(:,:,(levidx-1)*numOrien + orienidx) = strutA{levidx,orienidx};
		U(:,:,(levidx-1)*numOrien + orienidx) = strutU{levidx,orienidx};
		V(:,:,(levidx-1)*numOrien + orienidx) = strutV{levidx,orienidx};
		P(:,:,(levidx-1)*numOrien + orienidx) = strutP{levidx,orienidx};
		R(:,:,(levidx-1)*numOrien + orienidx) = ...
			sqrt(strutU{levidx,orienidx}.^2 + strutV{levidx,orienidx}.^2);
	end
end


% DCA: select pixel-wise whose A is largest
costFunc = A;
if(dca_opts.useWeightedAMFM == 1)
	for k=1:size(costFunc,3),
		costFunc(:,:,k) = costFunc(:,:,k) .* R(:,:,k);
	end
end

% use 3x3 window to smooth the cost function
if(dca_opts.useNeighbor == 1)
	h = fspecial('average', [dca_opts.window dca_opts.window]);
	for k=1:size(costFunc,3),
		costFunc(:,:, k) = imfilter(costFunc(:,:,k), h, 'replicate');
	end
end


% Select the maximum of the costFunc for each pixel
domA = zeros(M,N);
[domA, idx] = max(costFunc, [], 3);


% Find the values for 
domU = zeros(M,N);
domV = zeros(M,N);
domP = zeros(M,N);
for m=1:M
	for n=1:N
		domU(m,n) = U(m,n, idx(m,n));
		domV(m,n) = V(m,n, idx(m,n));
		domP(m,n) = P(m,n, idx(m,n));
	end
end
