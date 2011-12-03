% Perform Dominant component analysis from AM-FM components
% Perform pixel-wise dominant component with maximizing criteria A*R
%
% function [domA domU domV domP] = DCAFromComponents(A, U, V, P)
% 
%	inputs:
%			A: structure cell ( numLevels x numOrien x size(im) )
%			U: structure cell ( numLevels x numOrien x size(im) )
%			V: structure cell ( numLevels x numOrien x size(im) )
%			P: structure cell ( numLevels x numOrien x size(im) )
%
function [domA domU domV domP idx] = DCAFromComponents(...
									strutA, strutU, strutV, strutP, varargin);

if( nargin < 4 )
	error('Not enough input argument: DCA2(A, U, V, P)');
end

[numLevels, numOrien] = size(strutA);

% default values
useNeighbor = 0; % use 3x3 neighbor to select dominant amplitude
weightedAMFM = 0; % weighted AM with FM
numLevelsForDCA = numLevels; % just compute from 1-->K

dca_opts = [0 0 0];
if( length(varargin) == 1 )
	dca_opts = varargin{1};
end
useNeighbor = dca_opts(1);
weightedAMFM = dca_opts(2);
numLevelsForDCA = dca_opts(3);

if(numLevelsForDCA <= 0 || numLevelsForDCA > numLevels)
	%warning('The number of computed levels are out of range. Reset to 0');
	numLevelsForDCA = numLevels;
end


% Regenerate these structures into 3-D array
[M,N] = size(strutA{1,1});
A = zeros(M,N);
U = zeros(M,N);
V = zeros(M,N);
P = zeros(M,N);
R = zeros(M,N);
for levidx=1:numLevelsForDCA,
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
if(weightedAMFM == 1)
	for k=1:size(costFunc,3),
		costFunc(:,:,k) = costFunc(:,:,k) .* R(:,:,k);
	end
end

% use 3x3 window to smooth the cost function
if(useNeighbor == 1)
	h = fspecial('average', [11 11]);
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


%if(useNeighbor == 0) % use single pixel
%	for m=1:M,
%		for n=1:N,
%			tA = [];
%			for levidx=1:numLevelsForDCA,
%				for orienidx=1:numOrien,
%					if(weightedAMFM == 1)
%						tA = [tA (U{levidx,orienidx}(m,n)^2 + V{levidx,orienidx}(m,n)^2)*A{levidx,orienidx}(m,n)];
%					else
%						tA = [tA A{levidx,orienidx}(m,n)];
%					end
%				end
%			end
%
%			[y,iidx] = max(tA);
%			tlev = floor((iidx - 1) / numOrien) + 1;
%			torien = iidx - (tlev-1)*numOrien;
%
%			domA(m,n) = A{tlev, torien}(m,n);
%			domU(m,n) = U{tlev, torien}(m,n);
%			domV(m,n) = V{tlev, torien}(m,n);
%			domP(m,n) = P{tlev, torien}(m,n);
%		end
%	end
%else % use 3x3 neighbor for smoother results
%	for m=2:M-1,
%		for n=2:N-1,
%			tA = [];
%			% 3x3 neighborhood
%			for levidx=1:numLevelsForDCA,
%				for orienidx=1:numOrien,
%					if(weightedAMFM == 1)
%						aa = ((U{levidx,orienidx}(m-1:m+1,n-1:n+1)).^2 + ...
%							(V{levidx,orienidx}(m-1:m+1,n-1:n+1)).^2) .* A{levidx,orienidx}(m-1:m+1,n-1:n+1);
%						tA = [tA aa(:)'];
%					else
%						aa = A{levidx,orienidx}(m-1:m+1,n-1:n+1);
%						tA = [tA aa(:)'];
%						%tA = [tA A{levidx,orienidx}(m-1,n-1) A{levidx,orienidx}(m-1,n) A{levidx,orienidx}(m-1,n+1) ...
%								%A{levidx,orienidx}(m,n-1) A{levidx,orienidx}(m,n) A{levidx,orienidx}(m,n+1) ...
%								% A{levidx,orienidx}(m+1,n-1) A{levidx,orienidx}(m+1,n) A{levidx,orienidx}(m+1,n+1)]; 
%					end
%				end
%			end
%
%			[y,iidx] = max(tA);
%			tlev = floor((iidx - 1) / (9*numOrien)) + 1;
%			torien_a = iidx - (tlev-1)*(9*numOrien);
%			torien = floor( (torien_a - 1) / 9 ) + 1;
%			myidx = torien_a - (torien-1)*9;
%
%			if(myidx == 1)
%				mm = m-1;
%				nn = n-1;
%			elseif(myidx == 2)
%				mm = m;
%				nn = n-1;
%			elseif(myidx == 3)
%				mm = m+1;
%				nn = n-1;
%			elseif(myidx == 4)
%				mm = m-1;
%				nn = n;
%			elseif(myidx == 5)
%				mm = m;
%				nn = n;
%			elseif(myidx == 6)
%				mm = m+1;
%				nn = n;
%			elseif(myidx == 7)
%				mm = m-1;
%				nn = n+1;
%			elseif(myidx == 8)
%				mm = m;
%				nn = n+1;
%			else
%				mm = m+1;
%				nn = n+1;
%			end
%
%			domA(m,n) = A{tlev, torien}(mm,nn);
%			domU(m,n) = U{tlev, torien}(mm,nn);
%			domV(m,n) = V{tlev, torien}(mm,nn);
%			domP(m,n) = P{tlev, torien}(mm,nn);
%		end
%	end
%end %end if
