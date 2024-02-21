function [xCent, yCent] = find_center(image)
%Compute pupil centroid from SH lenslet image
%Used as preprocessing for wavefront reconstruction
% HISTORY:
% 2-19-2024 Warren Foster

    image = double(image);
    imSize = size(image);
    x = [1:imSize(1)];
    y = [1:imSize(2)];
    [X,Y] = meshgrid(y,x);
    imX = image .* X;
    imY = image .* Y;
    xCent = sum(sum(imX)) / sum(sum(image));
    yCent = sum(sum(imY)) / sum(sum(image));
end