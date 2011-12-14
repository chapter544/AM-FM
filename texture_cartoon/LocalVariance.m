function out = LocalVariance(in,L)

[M,N] = size(in);

out = zeros(M,N);
for m=L+1:1:M-1-L,
	for n=L+1:1:N-1-L
		a = in(m-L:m+L,n-L:n+L);
		out(m,n) = var(a(:));
	end
end
