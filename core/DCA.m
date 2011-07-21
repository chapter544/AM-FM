function [domA domU domV domP] = DCA(im, numLevels, numOrien, varargin)

if( nargin < 3 )
	error('Not enough input argument: DCA(im, numLevels, numOrien)');
end


% default values
useNeighbor = 0; % use 3x3 neighbor to select dominant amplitude
weightedAMFM = 0; % weighted AM with FM
excludedLevels = 0; % exclude from computing DCA
smoothingDemod = 0;

if( length(varargin) == 1 )
	useNeigbor = varargin{1};
end

if( length(varargin) == 2 )
	useNeigbor = varargin{1};
	excludedLevels = varargin{2};
end

if( length(varargin) == 3 )
	useNeigbor = varargin{1};
	excludedLevels = varargin{2};
	weightedAMFM = varargin{3};
end

if( length(varargin) == 4 )
	useNeigbor = varargin{1};
	excludedLevels = varargin{2};
	weightedAMFM = varargin{3};
	smoothingDemod = varargin{4};
end

if(excludedLevels < 0 || excludedLevels >= numLevels)
	warning('The number of excluded levels are out of range. Reset to 0');
	excludedLevels = 0;
end
numLevelsForDCA = numLevels - excludedLevels;




[M N] = size(im);
im = im - mean(im(:)); % eliminate DC component

% Generate the Steerable Pyramid filters
%[bandFilterDft, LoDft] = GenerateSteerablePyramid(M, N, numLevels, numOrien);
[bandFilterDft] = GenerateSteerablePyramidNoDC(size(im,1), size(im,2), numLevels, numOrien);

chanResponseR = cell(numLevels, numOrien);
chanResponseI = cell(numLevels, numOrien);

% Perform filtering through steerable pyramid
imDft = fftshift(fft2(im));
out = zeros(size(im));
for levidx = 1:numLevels,
	for orienidx = 1:numOrien,
		% Get the real part
		outDft = imDft .* bandFilterDft{levidx,orienidx};		
		chanResponseR{levidx,orienidx} = real(ifft2(ifftshift(outDft)));

		% Construct the imaginary part with the partial Hilbert transform
		H = GenerateDirectionalHilbertFilter(size(im,1), size(im,2), orienidx, numOrien);
		outImagDft = -sqrt(-1) .* outDft .* H; 
		chanResponseI{levidx,orienidx} = real(ifft2(ifftshift(outImagDft)));
	end
end
%Lo = real(ifft2(ifftshift(LoDft.*imDft)));


%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% DEMODUALTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
A = cell(numLevels,numOrien);
U = cell(numLevels,numOrien);
V = cell(numLevels,numOrien);
P = cell(numLevels,numOrien);
%Pls = cell(numLevels,numOrien);
%R = cell(numLevels,numOrien);
%T = cell(numLevels,numOrien);
if(smoothingDemod ~= 0)
	smoothingDemod = 1;
end
unwrapOptions = [300 smoothingDemod 5 5 0.5];
for levidx=1:numLevels,
	for orienidx=1:numOrien,
		[A{levidx,orienidx}, ...
			U{levidx,orienidx}, V{levidx,orienidx}, P{levidx,orienidx}] = ...
			phaseUnwrap(chanResponseR{levidx,orienidx}, ... 
							chanResponseI{levidx,orienidx}, unwrapOptions);

		%R{levidx,orienidx} = sqrt(U{levidx,orienidx}.^2 + V{levidx,orienidx}.^2);	
		%T{levidx,orienidx} = atan2(V{levidx,orienidx}, U{levidx,orienidx});
	end
end


% DCA: select pixel-wise whose A is largest
domA = zeros(M,N);
domU = zeros(M,N);
domV = zeros(M,N);
domP = zeros(M,N);
if(useNeighbor == 0) % use single pixel
	for m=1:M,
		for n=1:N,
			tA = [];
			for levidx=1:numLevelsForDCA,
				for orienidx=1:numOrien,
					if(weightedAMFM == 1)
						tA = [tA (U{levidx,orienidx}(m,n)^2 + V{levidx,orienidx}(m,n)^2)*A{levidx,orienidx}(m,n)];
					else
						tA = [tA A{levidx,orienidx}(m,n)];
					end
				end
			end

			[y,iidx] = max(tA);
			tlev = floor((iidx - 1) / numOrien) + 1;
			torien = iidx - (tlev-1)*numOrien;

			domA(m,n) = A{tlev, torien}(m,n);
			domU(m,n) = U{tlev, torien}(m,n);
			domV(m,n) = V{tlev, torien}(m,n);
			domP(m,n) = P{tlev, torien}(m,n);
		end
	end
else % use 3x3 neighbor for smoother results
	for m=2:M-1,
		for n=2:N-1,
			tA = [];
			% 3x3 neighborhood
			for levidx=1:numLevelsForDCA,
				for orienidx=1:numOrien,
					if(weightedAMFM == 1)
						aa = ((U{levidx,orienidx}(m-1:m+1,n-1:n+1)).^2 + ...
							(V{levidx,orienidx}(m-1:m+1,n-1:n+1)).^2) .* A{levidx,orienidx}(m-1:m+1,n-1:n+1);
						tA = [tA aa(:)'];
					else
						aa = A{levidx,orienidx}(m-1:m+1,n-1:n+1);
						tA = [tA aa(:)'];
						%tA = [tA A{levidx,orienidx}(m-1,n-1) A{levidx,orienidx}(m-1,n) A{levidx,orienidx}(m-1,n+1) ...
								%A{levidx,orienidx}(m,n-1) A{levidx,orienidx}(m,n) A{levidx,orienidx}(m,n+1) ...
								% A{levidx,orienidx}(m+1,n-1) A{levidx,orienidx}(m+1,n) A{levidx,orienidx}(m+1,n+1)]; 
					end
				end
			end

			[y,iidx] = max(tA);
			tlev = floor((iidx - 1) / (9*numOrien)) + 1;
			torien_a = iidx - (tlev-1)*(9*numOrien);
			torien = floor( (torien_a - 1) / 9 ) + 1;
			myidx = torien_a - (torien-1)*9;

			if(myidx == 1)
				mm = m-1;
				nn = n-1;
			elseif(myidx == 2)
				mm = m;
				nn = n-1;
			elseif(myidx == 3)
				mm = m+1;
				nn = n-1;
			elseif(myidx == 4)
				mm = m-1;
				nn = n;
			elseif(myidx == 5)
				mm = m;
				nn = n;
			elseif(myidx == 6)
				mm = m+1;
				nn = n;
			elseif(myidx == 7)
				mm = m-1;
				nn = n+1;
			elseif(myidx == 8)
				mm = m;
				nn = n+1;
			else
				mm = m+1;
				nn = n+1;
			end

			domA(m,n) = A{tlev, torien}(mm,nn);
			domU(m,n) = U{tlev, torien}(mm,nn);
			domV(m,n) = V{tlev, torien}(mm,nn);
			domP(m,n) = P{tlev, torien}(mm,nn);
		end
	end
end %end if
