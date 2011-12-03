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
trans_opt.withDCRing = 0;
trans_opt.phaseUnwrap = 1;
trans_opt.PHASE_CONST = 300;
trans_opt.smoothDemodulation = 1;
trans_opt.smoothWindow = 5;
trans_opt.smoothSigma = 0.1;
trans_opt.printParams = 0;
[A U V P Pls Resi] = AMFM_Transform(im_org, numLevels, numOrien, trans_opt);

anal_elapsed = toc(tanalysis); % elapsed time for analysis
fprintf('*** Analysis time: %0.2f seconds.\n\n', anal_elapsed);

tsynthesis = tic;
% reconstruct the image from AM-FM component
fprintf('-Perform reconstruction from the AM-FM model\n');
recon = Resi;
for levidx=1:numLevels,
	for orienidx=1:numOrien
		Precon = phaseReconLS(U{levidx,orienidx}, V{levidx,orienidx}, P{levidx,orienidx}(1,1), trans_opt.PHASE_CONST);
		recon = recon + A{levidx,orienidx} .* cos(Precon);
	end
end
syn_elapsed = toc(tsynthesis);
fprintf('*** Synthesis time: %0.2f seconds.\n\n', syn_elapsed);

% reconstruction error
[psnr, mse] = computePSNR(recon, im_org);
fprintf('MSE: %0.10f\n', mse);
fprintf('PSNR: %0.2f\n', psnr);


% computing the dominant component using DCA
dca_opts.useNeighbor = 0;
dca_opts.window = 5;
dca_opts.useWeightedAMFM = 0;
dca_opts.K = numLevels;
[dA dU dV dP dix] = DCAFromComponents(A, U, V, P, dca_opts);
