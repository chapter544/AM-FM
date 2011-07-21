function [H] = GenerateDirectionalHilbertFilter(row, col, orienidx, numOrien)

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


order = numOrien - 1;
lutsize = 1024;
factor = 2^(order) * factorial(order) / sqrt(numOrien*factorial(2*(order)));
Xcosn = pi * [-(3*lutsize+1):(lutsize+1)]/lutsize;
Ycosn = factor .* cos(Xcosn).^(order);

% Construct band-pass mask
k = orienidx;
interpOrg = Xcosn(1) + pi*(k-1)/numOrien + pi/(2*numOrien) + pi/2;
G = pointOp(angle, Ycosn, interpOrg, Xcosn(2) - Xcosn(1));

% Directional Hilbert Filter
H = -sign(G);

