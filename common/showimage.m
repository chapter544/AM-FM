function showimage(x, varargin)

inTitle = '';
if( size(varargin, 2) >= 1 )
	inTitle = varargin{1};
end

figure;
imagesc(x);
axis('image');
axis off
colormap(gray(256));
title(inTitle);

