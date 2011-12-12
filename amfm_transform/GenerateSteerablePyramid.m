%Function description:
%
%[bandFilterDft, LoDft] = GenerateSteerablePyramid(row, col, ... 
%													numLevels, numOrien)
%
%INPUTS:
%	row: row of image
%	col: column of image
%	numLevels: number of levels (decomposition)
%	numOrien: number of orientation per level
%
%OUTPUTS:
%	bandFilterDft: DFT of the Steerable Pyramid channel
% 
%Example: perform SP decomposition with 5 levels and 8 orientations
%	[bandFilterDft] = GenerateSteerablePyramid(im, 5, 8)
%
% Chuong Nguyen
% 07/07/09
% Modified 03/01/2010  to extends to rectangular images
% Modified 03/19/2010  changed to GenerateSteerablePyramid
%
%Credit: This implementation is based on the original 
%Simoncelli's steerable pyramid.
%			http://www.cns.nyu.edu/~eero/STEERPYR/
%

function [bandFilterDft, LoDft] = GenerateSteerablePyramid(row, col, numLevels, numOrien)

if nargin ~= 4 
	error('Not enough input arguments');
end


% Allocate memory  for outputs
bandFilterDft = cell(numLevels,numOrien);	% filter channel DFT
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
% Orientation part of the band-pass filter
% B = H * G where G is is the angular mask
order = numOrien - 1;
lutsize = 1024;
factor = 2^(order) * factorial(order) / sqrt(numOrien*factorial(2*(order)));
%Xcosn = pi * [-(2*lutsize+1):(lutsize+1)]/lutsize;
Xcosn = pi * [-(3*lutsize+1):(lutsize+1)]/lutsize;
Ycosn = factor .* cos(Xcosn).^(order);
for k=1:numOrien,
	% Construct band-pass filter
	interpOrg = Xcosn(1) + pi*(k-1)/numOrien + pi/(2*numOrien) + pi/2;
	G = pointOp(angle, Ycosn, interpOrg, Xcosn(2) - Xcosn(1), 0);
	bandFilterDft{1,k} = (-sqrt(-1))^(order) .* G .* hi0mask;
end


% recursively put low pass channel into the steerable pyramid
lo0mask = pointOp(log_rad, YIrcos, Xrcos(1), Xrcos(2)-Xrcos(1), 0);
loFilterDft{1} = lo0mask;

for levidx = 2:numLevels,
	Xrcos = Xrcos - log2(2);

	hiFilter = pointOp(log_rad, Yrcos, Xrcos(1), Xrcos(2)-Xrcos(1), 0);
	% Orientation part of the band-pass filter
	% B = H * G where G is is the angular mask
	for k=1:numOrien,
		% Construct band-pass filter
		interpOrg = Xcosn(1) + pi*(k-1)/numOrien + pi/(2*numOrien) + pi/2;
		G = pointOp(angle, Ycosn, interpOrg, Xcosn(2) - Xcosn(1), 0);
		bandFilterDft{levidx,k} = (-sqrt(-1))^(order) .* G .* ... .* 
									hiFilter .* loFilterDft{levidx-1};
	end

	% low-pass the residual image
	loFilter = pointOp(log_rad, YIrcos, Xrcos(1), Xrcos(2)-Xrcos(1), 0);
	loFilterDft{levidx} = loFilter .* loFilterDft{levidx-1};
end



% DO NOT break the low-pass into multiple orientations, since the 
% lowpass signal demodulation is error prone.
%for k=1:numOrien,
%	% Construct band-pass filter
%	interpOrg = Xcosn(1) + pi*(k-1)/numOrien + pi/(2*numOrien) + pi/2;
%	G = pointOp(angle, Ycosn, interpOrg, Xcosn(2) - Xcosn(1));
%	bandFilterDft{numLevels,k} = 
%				... (-sqrt(-1))^(order) .* G .* loFilterDft{numLevels-1}; 
%
%	% Perform filtering	
%	chanOutputDft{numLevels,k} = bandFilterDft{numLevels,k} .* imdft;
%end
LoDft = loFilterDft{numLevels}.^2;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SYNTHESIS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% put analysis channel output through complex conjugate filters (tight-frame)
for levidx = numLevels:-1:1,
	for orienidx = 1:numOrien,
		bandFilterDft{levidx,orienidx} = bandFilterDft{levidx,orienidx}.*conj(bandFilterDft{levidx,orienidx});
	end
end
