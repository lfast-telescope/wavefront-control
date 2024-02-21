function [arrows, idealCoords, actualCoords] = GetQuiver(image, referenceX, referenceY, magnification, output_plots)

% Get the quiver for actual points according to the regular grid.
% INPUTS:
%   image: captured image of the actual mirror;
%   referenceX/referenceY: x/y coordinate of the reference point;
%   magnification: magnification of the regular grid.
% OUTPUTS: 
%   arrows: N x 2 matrix showing displacements on the image;
%   idealCoords/actualCoords: N x 2 matrix showing coordinates of available
% ideal/actual points.
% HISTORY:
%   2023-11-02 - Yiyang Huang - initial implementation
%   2023-11-07 - Yiyang Huang - adjust the regular grid

% Generate the regular grid.
sizeX = size(image,2); sizeY = size(image,1);
xLeft = floor(referenceX/magnification);
xRight = floor((sizeX-referenceX)/magnification);
yLeft = floor(referenceY/magnification);
yRight = floor((sizeY-referenceY)/magnification);
xNum = xLeft + xRight + 1; yNum = yLeft + yRight + 1; % number of knots in x/y direction
x = magnification * (1:xNum) - magnification + 1;
y = magnification * (1:yNum) - magnification + 1;
displacementX = referenceX - x(xLeft+1);
displacementY = referenceY - y(yLeft+1);
x = x + displacementX; y = y + displacementY; % align with the reference point
[X,Y] = meshgrid(x,y); X = X(:); Y = Y(:); 

[xCent,yCent] = find_center(image);
radius = find_radius(image)*1.05; %Add a little fudge so that all the spots fit

if output_plots
    figure, imshow(image,[]);
    hold on
    scatter(X,Y,'r+');
    hold off
    exportgraphics(gca,'regular grid.png');
end

% Get the arrow for every available cross of the regular grid.
% Stretch contrast of the image.
image = TargetContrastControl(image); % improve uniformity 
image = double(image);
maxValue = max(image(:)); minValue = min(image(:));
image = uint8(255*(image-minValue)/(maxValue-minValue));
% Find arrows.
crossNum = length(X);
grayscaleThreshold = 24; % (Sometimes arbitrary numbers can eat your lunch)
halfSize = round(0.5*magnification); % (!)
availablePtsNum = 0; idealCoords = []; actualCoords = [];

if output_plots
    figure,imagesc(image)
    viscircles([xCent,yCent],radius)
end

%Evaluate if the cross should be evaluated within the pupil
for i = 1:crossNum
    % Get the to-be-judged region.
    lowerX = round(X(i)-halfSize); if lowerX < 1, lowerX = 1; end
    upperX = round(X(i)+halfSize); if upperX > sizeX, upperX = sizeX; end
    lowerY = round(Y(i)-halfSize); if lowerY < 1, lowerY = 1; end
    upperY = round(Y(i)+halfSize); if upperY > sizeY, upperY = sizeY; end
    targetRegion = image(lowerY:upperY,lowerX:upperX);
    % Judge if there is an available point.
    targetSizeX = size(targetRegion,2); targetSizeY = size(targetRegion,1);
    xTarget = 1:targetSizeX; yTarget = 1:targetSizeY; 
    [xTarget1D,yTarget1D] = meshgrid(xTarget,yTarget); 
    distance_from_center = sqrt((xTarget1D+lowerX-xCent).^2 + (yTarget1D+lowerY-yCent).^2);
    pixel_inside_circle = distance_from_center < radius;
    %If more than 80% of the region is outside the defined pupil, reject
    %this region
    if nnz(pixel_inside_circle) < numel(pixel_inside_circle) * 0.2
        
        if output_plots
            rectangle('Position',[lowerX lowerY (upperX-lowerX) (upperY-lowerY)],'EdgeColor','r')
        end
            continue
    else    
        xTarget1D = xTarget1D(:); yTarget1D = yTarget1D(:); % convert data to be 1D to faciliate judgment
        targetRegion1D = targetRegion(:);
        ptIndex = find(targetRegion1D==max(targetRegion1D));
        maxGrayscale = targetRegion1D(ptIndex);
        maxGrayscale = maxGrayscale(1); % ignore the impact of multiple maximum points
        if maxGrayscale > grayscaleThreshold % judge if there is a captured point
            xLocal = xTarget1D(ptIndex); yLocal = yTarget1D(ptIndex);
            xLocal = mean(xLocal); yLocal = mean(yLocal); % merge coordinates of multiple maximum points
            
            if output_plots
                rectangle('Position',[lowerX lowerY (upperX-lowerX) (upperY-lowerY)],'EdgeColor','c')
            end
        else
            if output_plots
                rectangle('Position',[lowerX lowerY (upperX-lowerX) (upperY-lowerY)],'EdgeColor','y')
            end
            continue;
        end
    end
    % Save the available point.
    xGlobal = xLocal + lowerX - 1; % go back into the original image
    yGlobal = yLocal + lowerY - 1;
    availablePtsNum = availablePtsNum + 1;
    idealCoords = [idealCoords; X(i) Y(i)];
    actualCoords = [actualCoords; xGlobal yGlobal];
end
arrows = actualCoords - idealCoords;

if output_plots
    figure, imshow(image); hold on
    quiver(idealCoords(:,1),idealCoords(:,2),arrows(:,1),arrows(:,2),'r');
    hold off
    exportgraphics(gca,'quiver.png');
end

end