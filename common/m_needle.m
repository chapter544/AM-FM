function [] = m_needle(u,v,weights,size,scale,proportion,density,ovlayBG)
% M_NEEDLE(U,V,WEIGHTES,SIZE,SCALE,PROPORTION,DENSITY,OVLAYBG)
%
% EXAMPLE USAGE:
%   m_needle(u,v,[],256,-1,'',8,[]);        Unit needle lengths
%   m_needle(u,v,[],256,0,'freq',8,[]);     Auto-scaled frequency
%   m_needle(u,v,[],256,10,'freq',8,[]);    User-scaled frequency
%   m_needle(u,v,a,256,0,'freq',8,[]);      Auto-scaled weighted frequency
%   m_needle(u,v,a,256,0.1,'freq',8,[]);    User-scaled weighted frequency
%   m_needle(u,v,[],256,-1,'',8,a);         Unit needle overlay

% extract the needed frequency values
u = u(density:density:size-1,density:density:size-1);
v = v(density:density:size-1,density:density:size-1);

% compute the needle angles
theta = atan2(v,u);

% compute the needle magnitudes
if(scale < 0)
    % make all needle magnitudes equal to the needle spacing
    needleMag = density;
else
    % compute the frequency magnitudes
    needleMag = sqrt(u.^2 + v.^2);
    % make the needle magnitudes proportional to wavelength
    if(strcmp(upper(proportion),'WAVE'))
        smallMags = find(needleMag < 0.01);
        needleMag(smallMags) = 1;
        needleMag = 1./needleMag;
        needleMag(smallMags) = 0;
    elseif(strcmp(upper(proportion),'FREQ'))
        % nothing to do here
    else
        error('Unrecognized proportion type ''%s''',proportion);
    end
    % weight the needle magnitudes when weights are specified
    if(~isempty(weights))
        weights = weights(density:density:size-1,density:density:size-1);
        needleMag = needleMag.*abs(weights);
    end
    % scale the needle magnitudes (automatically or by specified constant)
    if(scale == 0)
        % auto-scale
        medMag = medfilt2(needleMag);
        avgMedMag = mean(mean(medMag));
        stdMedMag = sqrt(mean(mean((medMag - avgMedMag).^2)));
        medMagIndices = find(abs(medMag - avgMedMag) < 3.*stdMedMag);
        maxMag = max(max(medMag(medMagIndices)));
        needleMag = needleMag.*density./maxMag;
    else
        needleMag = needleMag.*scale;
    end
end

% compute the scaled frquencies
u = needleMag.*cos(theta);
v = needleMag.*sin(theta);

figure;
% if making an overlay, plot the background
if(~isempty(ovlayBG))
    % stretch background to gray scale range
    minBG = min(min(ovlayBG));
    maxBG = max(max(ovlayBG));
    ovlayBG = 240.*(ovlayBG - minBG);
    if(maxBG ~= minBG)
        ovlayBG = ovlayBG./(maxBG - minBG);
    end
    ovlayBG = round(ovlayBG);
    image(ovlayBG);
    colormap(gray(256));
    % overlay the needle plot
    hold on;
    [x,y] = meshgrid(density:density:size-1,density:density:size-1);
    quiver(x,y,u,v,0,'-w');
    hold off;
    % fix the figure orientation, size, aspect ratio, and border
    axis('ij'); axis([1,size,1,size]); axis('square'); axis('off');
    fixFigure(gcf,1.5.*size,0,0,0,0);
else
    % make the needle plot
    [x,y] = meshgrid(density:density:size-1,density:density:size-1);
    quiver(x,y,u,v,0,'-k');
    % fix the figure orientation, size, aspect ratio, and border
    axis('ij'); axis([1,size,1,size]); axis('square'); axis('off');
    fixFigure(gcf,1.5.*size,5,5,5,5);
end
