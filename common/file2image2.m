function [img] = file2image2(dataType,row,col,fileName)
%FILE2IMAGE Read a raw image from disk
%  [IMG] = FILE2IMAGE(DATATYPE,row,col,FILENAME) returns the SIZExSIZE
%  raw image specified by FILENAME and DATATYPE. See the FREAD help
%  section for common values of DATATYPE. Here are some examples:
%
%  DATATYPE     Bits Per Pixel
%  --------     --------------
%  'uchar'      8  (unsigned)
%  'float'      32
%  'double'     64


% read in the image
fileId = fopen(fileName,'r');
[img numPix] = fread(fileId,[row,col],dataType);
fclose(fileId);

% check for formatting error
if(numPix ~= row*col)
    disp('Error reading image');
    return;
end

% fix Matlab transpose of image
img = img';

return;
