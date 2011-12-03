clear all
close all
clc

M = 256;
N = 256;

%im_org = file2image2('float', M, N, sprintf('Images/%sR.float', rootName));
im_org = generateChirp(M, N);

tanalysis = tic;

numLevels = 5;
numOrien = 8;
trans_opt = [0 1 300 1 3 0.5];
[A U V P Pls Resi] = AMFM_Transform(im_org, numLevels, numOrien);

anal_elapsed = toc(tanalysis); % elapsed time for analysis
fprintf('*** Analysis time: %0.2f seconds.\n\n', anal_elapsed);

tsynthesis = tic;
% reconstruct the image from AM-FM component
fprintf('-Perform reconstruction from the AM-FM model\n');
recon = Resi;
for levidx=1:numLevels,
	for orienidx=1:numOrien
		Precon = phaseReconLS(U{levidx,orienidx}, V{levidx,orienidx}, P{levidx,orienidx}(1,1), trans_opt(3));
		recon = recon + A{levidx,orienidx} .* cos(Precon);
	end
end
syn_elapsed = toc(tsynthesis);
fprintf('*** Synthesis time: %0.2f seconds.\n\n', syn_elapsed);

% reconstruction error
[psnr, mse] = computePSNR(recon, im_org);
fprintf('MSE: %0.10f\n', mse);
fprintf('PSNR: %0.2f\n', psnr);
