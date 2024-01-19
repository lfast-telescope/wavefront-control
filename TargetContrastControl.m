function correctedIm = TargetContrastControl(image)

% This function is used to pre-process the target images. Primary work is 
% to improve image contrast.
% HISTORY:
%   2023-08-08 - Yiyang Huang - initial implementation

% Stretch the grayscale.
image = double(image);
maxValue = max(image(:)); minValue = min(image(:));
transIm = uint8(255*(image-minValue)/(maxValue-minValue));
% figure, imshow(transIm);
% pause(0.2);

% Improve imaging effects.
correctedIm = adapthisteq(transIm); % adaptive local histograms
% figure, imshow(correctedIm);
% pause(0.2);
% close all;

end
