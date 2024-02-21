function [radius] = find_radius(image)

%Compute pupil radius from SH lenslet image
%Used as preprocessing for wavefront reconstruction
% HISTORY:
% 2-19-2024 Warren Foster


    image = double(image);
    maxValue = max(image(:)); minValue = min(image(:));
    imageUint8 = uint8(255*(image-minValue)/(maxValue-minValue));
    binaryIm = imbinarize(imageUint8,'global'); % global thresholding
    [row,col] = find(binaryIm);
    rowD = max(row) - min(row);
    colD = max(col) - min(col);
    radius = (rowD + colD) / 4;
end