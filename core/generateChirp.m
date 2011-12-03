function [f] = generateChirp(M,N)
f = zeros(M,N);
for m=1:M
	for n=1:N
		%f(m,n) = cos( 12^2*((m-M/2)^2 + (n-N/2)^2) / (M*N) );
		%a(m,n) = exp( - ((m-M/2)^2 + (n-N/2)^2) / (M*N) );
		a = exp( - ((m-M/2)^2 + (n-N/2)^2) / (M*N) );
		p = cos( 12^2*((m-M/2)^2 + (n-N/2)^2) / (M*N) );
		f(m,n) = a * p;
	end
end
f = f - mean(mean(f));
%file2image(f, 'float', 'chirpR.float');
