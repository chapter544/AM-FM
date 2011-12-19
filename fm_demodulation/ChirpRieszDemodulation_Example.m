%function [texture, cartoon] = MyPu4(f, a_threshold, b_threshold)

clear all
close all
clc

%M = 512;
%N = 512;
M = 256;
N = 256;
%rootName = 'barbara512R';
%rootName = 'chirpR';
%rootName = 'lenaR';
%rootName = 'burlapR';
%rootName = 'mandrilR';
%rootName = 'barbara_face256R';

%f  = file2image2('float', M, N, sprintf('Images/%s.float', rootName));

f = zeros(M,N);
a_org = zeros(M,N);
p_org = zeros(M,N);
alpha = 6;
beta = 8;
for m=1:M
	for n=1:N
		a_org(m,n) = exp( - 2*((m-M/2)^2 + (n-N/2)^2) / (M*N) );
		p_org(m,n) = (alpha * 2*pi*(m - M/2)^2 / M^2 + beta * 2*pi*(n-N/2)^2/ N^2) * 2;
	end
end
f = a_org .* cos(p_org);
f = f - mean(mean(f));

%showimage(a_org);
%showimage(cos(p_org));
%showimage(f);


[a, u, v, p, pls] = RieszPhaseUnwrap(f, 300.0, 1);

%r2 = sqrt(u2.^2 + v2.^2);
%figure;
%hist(r2(:), 100);

%showimage(a2, 'a2');
%showimage(cos(p2), 'cos(p2)');
close all;

L = 1:1:N;
row = 128;
figure;
plot(L,f(row,:), 'r'); title('original');
hold on
plot(L,a(row,:), 'b'); title('demodulated a');
hold on
plot(L,cos(p(row,:)), 'c'); title('demodulated cosP');
hold off
legend('original', 'a', 'cosP');
axis([1 N -1.2 1.2]);



figure;
subplot(2,3,1)
imagesc(f); axis('image'); axis off; colormap(gray(256)); title('Original');
subplot(2,3,2)
imagesc(a_org); axis('image'); axis off; colormap(gray(256)); title('Aorg');
subplot(2,3,3)
imagesc(cos(p_org)); axis('image'); axis off; colormap(gray(256)); title('cos(Porg)');
subplot(2,3,4)
%imagesc(a_org - a2); axis('image'); axis off; colormap(gray(256)); title('a_org - a diff');
mesh(a_org - a); axis('image'); axis off; colormap(gray(256)); title('a_org - a diff');
subplot(2,3,5)
imagesc(a); axis('image'); axis off; colormap(gray(256)); title('a demod');
subplot(2,3,6)
imagesc(cos(p)); axis('image'); axis off; colormap(gray(256)); title('cosP demod');



%
%t = -pi/2:pi/100:pi/2;
%y = 2*( 1 - 1 ./ (1 + exp( - 1* abs(t - pi/4) )));
%plot(t, y);
%
%
%
%
%close all;
%alpha = 0.6;
%beta = 0.8;
%gamma_r = 0.01;
%gamma_t = 1;
%t_filt = pi/4;
%
%texture = zeros(size(f));
%for levidx=1:numLevels
%	%gain = 1 ./ (1 + exp( - gamma * (R{levidx} - (alpha+beta)/2) ) );
%	gain_t = 2 * ( 1 - 1 ./ (1 + exp( - gamma_t * abs(abs(T{levidx}) - t_filt))) );
%	texture = texture + gain_t .* A{levidx} .* cos(P{levidx});
%end
%subplot(1,2,1)
%imagesc(texture); axis('image'); axis off; colormap(gray(256));
%subplot(1,2,2)
%imagesc(f-texture); axis('image'); axis off; colormap(gray(256));
%















