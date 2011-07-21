function [bandFilterDft] = Generate1LevelSteerableFilters(row, col, numOrien)

if nargin ~= 3 
	error('Not enough input arguments');
end


% Allocate memory  for outputs
bandFilterDft = cell(1,numOrien);	% filter channel DFT

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
	bandFilterDft{1,k} = (-sqrt(-1))^(order) .* G;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SYNTHESIS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% put analysis channel output through complex conjugate filters (tight-frame)
for orienidx = 1:numOrien,
	bandFilterDft{1,orienidx} = bandFilterDft{1,orienidx}.*conj(bandFilterDft{1,orienidx});
end
