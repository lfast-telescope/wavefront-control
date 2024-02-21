clearvars; clc; close all;

file_name = 'C:\Users\warrenbfoster\OneDrive - University of Arizona\Documents\LFAST\on-sky\20240118\20240118T0228_mirfak_tecsoff.bmp';
output_plots = false;

%function [highOrderSurface] = wavefrontReconstruction(im_filename,output_plots, zernikeTerms)

% Read the image.
imageColor = imread(file_name);
image = rgb2gray(imageColor);

%Crop image: locate pupil and center on array
imCrop = CropImage(image);

if output_plots
    figure, imshow(image,[]);
end

%Determine rotation correction to make spot field align with sensor array 
angle_start = 0;
angle_range = 1;

for trial = [1:5]
     best_angle = find_best_rotation(imCrop,angle_start,angle_range);
     if (best_angle == angle_start + angle_range) || (best_angle == angle_start - angle_range)
     else
         angle_range = angle_range / 2;
     end
     angle_start = best_angle;
end

image = imrotate(imCrop,best_angle);


% Get parameters of the regular grid according to known bright spots.
[referenceX, referenceY, magnification] = GetGrid(image, output_plots);

if false
    magnification = 54; % magnification should be a constant for the same Shack-Hartmann wavefront sensor
end
% Calculate the quiver.
[arrows, idealCoords, ~] = GetQuiver(image, referenceX, referenceY, magnification, output_plots);

% Calibrate slope data according to the quiver.
s = 5.63116; % distance between 2 stars in pixels
tanTheta = 5.5; % angular separation of 2 stars in arcseconds
tanTheta = tanTheta * 4.848e-6; % angular separation in radians
d = s / tanTheta; % calibrated distance in pixels between the lens array and sensor
slope = 0.5 * arrows / d; % calibrated slope

% Get the regular slope maps.
ptsNumX = (max(idealCoords(:,1)) - min(idealCoords(:,1))) / magnification;
ptsNumY = (max(idealCoords(:,2)) - min(idealCoords(:,2))) / magnification;
N = (ptsNumX + ptsNumY) / 2; % number of separations of micro lenses on the diameter
slopeMagnification = 1; % equal to 1 once the slope is calibrated in previous sections
mirrorDiameter = 30 * 25.4; % diameter of the primary mirror in millimeters
actualSpacing = mirrorDiameter / N; % actual spacing on the primary mirror indicated by micro lenses
lateralMagnification = actualSpacing / magnification; % actual spacing on the primary mirror covered by a pixel

%Start here tomorrow

integrationStep = 3;
[regularSlopeX,regularSlopeY,xCoordinates,yCoordinates] = ...
    Quiver2RegularSlope(slope,idealCoords,slopeMagnification, ...
    lateralMagnification,integrationStep);
% yCoordinates = flip(yCoordinates); % correct the flipped y axis of the pixel coordinate system

% Display and save the regular slope maps.
mirrorPos = lateralMagnification * idealCoords;
zCoords = 50 * ones(size(mirrorPos,1),1);

% Integrate slope maps.
shapeDiff = SlopeIntegration(regularSlopeX, regularSlopeY, 'Southwell');

if output_plots
    figure, mesh(xCoordinates,yCoordinates,regularSlopeX); hold on
    scatter3(mirrorPos(:,1),mirrorPos(:,2),zCoords,'k')
    colormap('jet'); colorbar; view([0 90]); axis equal
    title('Regular slope x');
    set(gca,'YDir','reverse'); % set the inverted y axis to recover original display effect
    hold off

    exportgraphics(gca,['regular slope X','.png']);
    figure, mesh(xCoordinates,yCoordinates,regularSlopeY); hold on
    scatter3(mirrorPos(:,1),mirrorPos(:,2),zCoords,'k')
    colormap('jet'); colorbar; view([0 90]); axis equal
    title('Regular slope y');
    set(gca,'YDir','reverse'); % set the inverted y axis to recover original display effect
    hold off
    exportgraphics(gca,['regular slope Y','.png']);

    figure, imhandle = meshc(xCoordinates,yCoordinates,shapeDiff); 
    imhandle(2).EdgeColor = 'k'; imhandle(2).ZLocation = 'zmax';
    imhandle(2).LevelStep = 250;
    colormap('jet'); colorbar; view([0 90]); axis equal
    title('Shape difference');
    % caxis([-2000 2000]);
    set(gca,'YDir','reverse'); % set the inverted y axis to recover original display effect
    exportgraphics(gca,['shape difference','.png']);
end
%%
% Remove low-order terms.
xCentered = xCoordinates - mean(xCoordinates(:));
yCentered = yCoordinates - mean(yCoordinates(:));
yCentered = flip(yCentered);
zernikeTerms = [1,2,3,4];
[highOrderMap,fitCoeff] = RemoveLowOrderZernike(xCentered, yCentered, shapeDiff, zernikeTerms, output_plots);
