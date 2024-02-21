function [image_crop] = CropImage(image)

% Get parameters of the regular grid.
% INPUT:
%   image: captured image of the actual mirror;
% OUTPUTS: 
%   referenceX/referenceY: x/y coordinate of the reference point;
%   magnification: magnification of the regular grid.
% HISTORY:
%   2-19-2024: Warren Foster

% Characterize SH pupil
image = double(image);
[xCent,yCent] = find_center(image);
xCent = uint16(xCent);
yCent = uint16(yCent);
radius = find_radius(image)*1.05; %Add a little fudge so that all the spots fit

if false
    imagesc(image)
    viscircles([xCent,yCent], radius,'Color','r');
end

xSize = size(image,2); ySize = size(image,1);
yTop = floor(yCent);
yBot = floor((ySize-yCent));
xLeft = floor(xCent);
xRight = floor((xSize-xCent));

if yTop < yBot
    image_crop = image(1:(yCent+yTop),:);
else
    image_crop = image((yCent-yBot):end,:);
end

if xLeft < xRight
    image_crop = image_crop(:,1:(xCent+xLeft));
else
    image_crop = image_crop(:,(xCent-xRight):end);
end

end
