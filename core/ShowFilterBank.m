function ShowFilterBank(M,N)

numLev = 5;
numOrien = 8; 
[bandFilterDft] = GenerateSteerablePyramidNoDC(M,N, numLev, numOrien);

out = zeros(M,N);
for idx=1:numLev,
	for idx2=1:numOrien,
		temp = abs(bandFilterDft{idx,idx2});
		for i=1:M,
			for j=1:N,
				if(out(i,j) < temp(i,j))
					%out(i,j) = out(i,j) + temp(i,j);	
					out(i,j) = temp(i,j);	
				end
			end
		end
	end
end

image2file(out, 'float', 'steerablePyramid.float', 0);


