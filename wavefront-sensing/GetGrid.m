function [referenceX, referenceY, magnification] = GetGrid(image,output_plots)

% Get parameters of the regular grid.
% INPUT:
%   image: captured image of the actual mirror;
% OUTPUTS: 
%   referenceX/referenceY: x/y coordinate of the reference point;
%   magnification: magnification of the regular grid.
% HISTORY:
%   2023-12-13 - Yiyang Huang - initial implementation
%   2024-01-16 - Warren Foster - add output plot suppression
%   2024-01-18 - Yiyang Huang - change the way of determining magnification


% Get pixel coordinates of bright spots.

% Binarize the image to get a mask.
image = double(image);
maxValue = max(image(:)); minValue = min(image(:));
imageUint8 = uint8(255*(image-minValue)/(maxValue-minValue));
binaryIm = imbinarize(imageUint8,'global'); % global thresholding

if output_plots
    figure, imshow(binaryIm,[]);
end

% Use the mask to get the target regions on the smoothed image.
sigma = 2; % (!)
smoothIm = imgaussfilt(image,sigma);
maskedIm = binaryIm .* smoothIm;

if output_plots
    figure,imshow(maskedIm,[]);
end

% Find local maximums of target regions to show positions of bright spots.
brightSpots = imregionalmax(maskedIm);

if output_plots
    figure, imshow(brightSpots);
    close all;
end

[rows,cols] = find(brightSpots == 1);

% Get magnification and central reference of bright spots.
% Calculate the magnification.
spotsNum = length(rows);
spotSet = [cols rows];
for i = 1:spotsNum
    judgedSpot = spotSet(i,:);
    judgedSpotSet = spotSet; judgedSpotSet(i,:) = [];
    [closestSpot, ~] = FindClosestPoint(judgedSpotSet, judgedSpot);
    distance(i) = norm(closestSpot-judgedSpot);
end
sortedDist = sort(distance);
croppedPercent = 0.1; % (!)
croppedNum = round(croppedPercent*spotsNum);
croppedDist = sortedDist(croppedNum+1:spotsNum-croppedNum);
magnification = max(croppedDist);
% Find the central reference.
rowAve = mean(rows); colAve = mean(cols);
centralSpot = [colAve rowAve];
[closestSpot, ~] = FindClosestPoint(spotSet, centralSpot);
referenceX = closestSpot(1);
referenceY = closestSpot(2);

if output_plots
    figure,imshow(maskedIm);
    hold on;
    scatter(centralSpot(1),centralSpot(2),'Color','g')
    scatter(closestSpot(1),closestSpot(2),'r')
    hold off

end
