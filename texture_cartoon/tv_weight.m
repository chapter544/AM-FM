function y = tv_weight(x)

[M,N] = size(x);
y = zeros(size(x));

a1 = 0.25;
a2 = 0.5;

for m=1:M,
	for n=1:N,
		if(x(m,n) < a1);
			y(m,n) = 0.0;
		elseif( x(m,n) > a2);
			y(m,n) = 1.0;
		else
			y(m,n) = 1 / (a2 - a1) * x(m,n) - 1;
		end
	end
end
