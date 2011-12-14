%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FullScaleStretch2.m
% [y] = FullScaleStretch(x, varargin)
% @ x       :   input 2-D image
% @ y       :   full-scale stretch of image
%
% Assume that the image is 8bit/pixel
% NTC 3/1/2006
% Modified 4/12/2006
%   Fix the non-integer stretch of unit8
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [y] = FullScaleStretch2(x, varargin)

% Obtain row and column of image
m_row = size(x,1);
m_col = size(x,2);

% Obtain maximum and minimum values image
if( length(varargin) == 0)
	m_imMax = max(max(x));
	m_imMin = min(min(x));
elseif( length(varargin) == 2) 
	m_imMax = varargin{1};
	m_imMin = varargin{2};
else 
	error('not enough arguments');
end

% Check to make sure m_imMax > m_imMin
if (m_imMax <= m_imMin)
    y = zeros(m_row,m_col);    
else
    % Scale to full dynamic
    m_scaleFactor = 255/(m_imMax - m_imMin);
    y = uint8(floor(m_scaleFactor .* (x - m_imMin) + 0.5));
end
