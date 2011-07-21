function [] = image2file(img,dataType,fileName,cs)
%IMG2FILE Write (contrast stretched) image to disk
%  [] = IMAGE2FILE(IMG,DATATYPE,FILENAME,CS) writes IMG to FILENAME
%  using DATATYPE for each pixel. If DATATYPE is set to 'jpg', then
%  the image is contrast stretched and written as JPEG. For all other
%  values of DATATYPE, the image is contrast stretched before writing
%  when the CS flag is nonzero. See the FREAD help section for common
%  values of DATATYPE. Here are some examples:
%
%  DATATYPE     Bits Per Pixel
%  --------     --------------
%  'uchar'      8  (unsigned)
%  'float'      32
%  'double'     64


if(strcmp(upper(dataType),'JPG'))
    % contrast stretch image to range (0,255)
    img = contStretch(img);
    % save the image as a jpeg
    imwrite(img + 1,gray(256),fileName,'jpg','quality',100);
else
    if(cs ~= 0)
        % contrast stretch image to range (0,255)
        img = contStretch(img);
    end
    % fix Matlab transpose of image
    img = img';
    % save the image in raw format
    fileId = fopen(fileName,'w');
    fwrite(fileId,img,dataType);
    fclose(fileId);
end

return;