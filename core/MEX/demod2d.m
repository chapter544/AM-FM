function [a, u, v, p] = demod2d(img, numRows, numCols)
%  Function demodulates a complex image into its AM and FM components. The FM
%    component is computed using a tensor product cubic spline model applied to
%    the scaled, unwrapped phase embedded with congruent principal values. This
%    unwrapped phase modulation is also returned. At pixels where the amplitude
%    is zero, the phase is copied from a neighboring pixel with a valid
%    principal value. Returns on error.
%
%  10/21/2005  ras

numPix = numRows*numCols;
numZeroAM = 0;
zeroAM = zeros(numRows,numCols);
validPixel = zeros(numRows,numCols);

%DBL_EPSILON = 0.00001;
PHASE_SCALE = 300;

%compute the amplitude and phase modulations
a = abs(img);
p = angle(img);

for row = 1:numRows
    for col = 1:numCols
        if (abs(real(img(row,col))) == 0) && (abs(imag(img(row,col))) == 0)
            %set zeroAM flag and increment its counter
            zeroAM(row,col) = 1;
            numZeroAM = numZeroAM + 1;
        else
            %clear zeroAM flag
            zeroAM(row,col) = 0;
            validPixel(row,col) = 1;

            %calculate wrapped phase in (-pi,pi]
            if p(row,col) == -pi
                p(row,col) = pi;
            end
        end
    end 
end

%return when amplitude is zero everywhere
if (numZeroAM == numPix)
    return;
end

% interpolate the phase at pixels with zero amplitude by copying a valid
%   phase value from a neighboring pixel; neighbors are searched in the order
%   left, right, top, bottom, and mirror symmetry is used at the boundaries
while (numZeroAM > 0)
    for row = 1:numRows
        for col = 1:numCols
            if ( zeroAM(row,col) ~= 0 )
                
                % process the left phase value
                if (col == 1) 
                    if(zeroAM(row,(col+1)) == 0) 
                        p(row,col) = p(row,(col+1));
                        zeroAM(row,col) = 0;
                        numZeroAM = numZeroAM - 1;
                        break;
                    end
                else
                    if(zeroAM(row,(col-1)) == 0) 
                        p(row,col) = p(row,(col-1));
                        zeroAM(row,col) = 0;
                        numZeroAM = numZeroAM - 1;
                        break;
                    end
                end
                % process the right phase value
                if (col == numCols) 
                    if(zeroAM(row,(col-1)) == 0) 
                        p(row,col) = p(row,(col-1));
                        zeroAM(row,col) = 0;
                        numZeroAM = numZeroAM - 1;
                        break;
                    end
                else
                    if(zeroAM(row,(col+1)) == 0) 
                        p(row,col) = p(row,(col+1));
                        zeroAM(row,col) = 0;
                        numZeroAM = numZeroAM - 1;
                        break;
                    end
                end
                % process the top phase value
                if (row == 1) 
                    if(zeroAM((row+1),col) == 0) 
                        p(row,col) = p((row+1),col);
                        zeroAM(row,col) = 0;
                        numZeroAM = numZeroAM - 1;
                        break;
                    end
                else
                    if(zeroAM((row-1),col) == 0) 
                        p(row,col) = p((row-1),col);
                        zeroAM(row,col) = 0;
                        numZeroAM = numZeroAM - 1;
                        break;
                    end
                end
                % process the bottom phase value
                if (row == numRows) 
                    if(zeroAM((row-1),col) == 0) 
                        p(row,col) = p((row-1),col);
                        zeroAM(row,col) = 0;
                        numZeroAM = numZeroAM - 1;
                        break;
                    end
                else
                    if(zeroAM((row+1),col) == 0) 
                        p(row,col) = p((row+1),col);
                        zeroAM(row,col) = 0;
                        numZeroAM = numZeroAM - 1;
                        break;
                    end
                end

            end
        end
    end
end


%unwrap the phase modulation
p = unwrap2d(p,numRows,numCols);

%compute the phase gradient using a tensor product cubic spline model
[u,v] = splineGrad(p,numRows,numCols,1,1);

%remove the phase scale from the gradient
u = u / PHASE_SCALE;
v = v / PHASE_SCALE;
