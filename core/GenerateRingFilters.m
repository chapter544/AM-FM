% Generate Ring Filters (like donought type)
%function [bandFilterDft] = GenerateRingFilters(row, col, numLevels)
%	
% inputs:
%			row: vertical dimension
%			col: horizontal dimension
%			numLevels: number of levels (scales)
%
% outputs:
%			bandFilterDft: cell structure DFT of the filter scales. 
%				For example, bandFilterDft{1} is the highest freq. level
%
%
function [bandFilterDft] = GenerateRingFilters(row, col, numLevels)

if nargin ~= 3 
	error('Not enough input arguments');
end



bandFilterDft = cell(numLevels,1);	% filter channel DFT
loFilterDft = cell(numLevels,1);			% low pass filter of each level

dims = [row col];
ctr = ceil((dims+0.5)/2);
[xramp,yramp] = meshgrid( ([1:dims(2)]-ctr(2))./(dims(2)/2), ... 
							([1:dims(1)]-ctr(1))./(dims(1)/2) );

% radial portion in log scale
log_rad = sqrt(xramp.^2 + yramp.^2);
log_rad(ctr(1),ctr(2)) =  log_rad(ctr(1),ctr(2)-1);
log_rad  = log2(log_rad);

% angle portion
angle = atan2(yramp,xramp);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Filter interpolation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sz = 256;
Xrcos = [-sz-1:1] / sz;
X_rad = 0.5 * pi .* Xrcos ; % X_rad = [-pi/2, 0]
Yrcos = cos(X_rad).^2;
% Symmetry for interpolation
Yrcos(1) = Yrcos(2);
Yrcos(sz+3) = Yrcos(sz+2);
Yrcos = sqrt(Yrcos); % sqrt of coinse^2
YIrcos = sqrt(1 - Yrcos.^2);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ANALYSIS 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% first hi-pass channel
hi0mask = pointOp(log_rad, Yrcos, Xrcos(1), Xrcos(2)-Xrcos(1), 0);
bandFilterDft{1} = hi0mask;

% recursively put low pass channel into the steerable pyramid
lo0mask = pointOp(log_rad, YIrcos, Xrcos(1), Xrcos(2)-Xrcos(1), 0);
loFilterDft{1} = lo0mask;

for levidx = 2:numLevels-1,
	Xrcos = Xrcos - log2(2);

	hiFilter = pointOp(log_rad, Yrcos, Xrcos(1), Xrcos(2)-Xrcos(1), 0);
	bandFilterDft{levidx} = hiFilter .* loFilterDft{levidx-1};

	% low-pass the residual image
	loFilter = pointOp(log_rad, YIrcos, Xrcos(1), Xrcos(2)-Xrcos(1), 0);
	loFilterDft{levidx} = loFilter .* loFilterDft{levidx-1};
end

% DC level
bandFilterDft{numLevels} = loFilterDft{numLevels-1};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SYNTHESIS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% put analysis channel output through complex conjugate filters (tight-frame)
for levidx = numLevels:-1:1,
	bandFilterDft{levidx} = bandFilterDft{levidx}.*conj(bandFilterDft{levidx});
end
