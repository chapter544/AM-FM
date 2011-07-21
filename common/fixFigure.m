function fixFigure(h,size,top,bot,left,right)
% FIXFIGURE(H,SIZE,TOP,BOT.LEFT,RIGHT) Resize figure and its border
%  Function makes figure with handle H a SIZExSIZE pixel figure and adds a
%  border of TOP, BOT, LEFT, and RIGHT pixels.

set(h,'units','pixels'); p = get(h,'position');
set(h,'position',[p(1) p(2) size+left+right size+top+bot]);
a = get(h,'currentaxes');
set(a,'units','pixels');
set(a,'position',[left+1 bot+1 size size]);