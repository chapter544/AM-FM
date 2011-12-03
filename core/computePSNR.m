% [psnr] = computePSNR(x, y)
% Compute the Peak SNR of two signal
% psnr = 10*log10( max(x)^2 / mse )
% x: original (good image)
% y: original (reconstructed image)
function [psnr, mse] = computePSNR(x,y)
% convert to double precision
xx = double(x);
yy = double(y);

%max_val = max( max(xx(:)), max(yy(:)) );
max_val = max( max(xx(:)) );
mse = mean( (xx(:) - yy(:)).^2 );
psnr = 10*log10(max_val^2 / mse);
