% Plot the needle diagram with U and V
% A (if provided) is used as a weight matrix
% for the FM function
function myNeedlePlot(A, U, V, L, gain)
% check for image dimension argreement
if( size(U,1)  ~= size(V,1) || ... 
	size(U,2)  ~= size(V,2) )
	error('Real and imaginary image are of different size');
end

if( size(A,1)  == 0 || size(A,2) == 0)
	A = ones(size(U));
	% fix the boundary
	A(1,:) = 0;
	A(end,:) = 0;
	A(:,1) = 0;
	A(:,end) = 0;
end

% fix the boundary
[M,N] = size(U);
figure;
imagesc(A); axis('image'); colormap(gray(256)); axis off;
hold on;
xxidx = L:L:N-L;
yyidx = L:L:M-L;
[xx,yy] = meshgrid(xxidx,yyidx);
AA = A(xxidx, yyidx);
UU = AA.*U(xxidx, yyidx);
VV = AA.*V(xxidx, yyidx);
quiver(xx,yy, UU, VV, gain,'color', 'r', 'linewidth', 1);
hold off;
