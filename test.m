clearvars; clc; close all;

file_name = 'C:\Users\warre\OneDrive\Documents\LFAST\on-sky\2024-01-10\WFS\20240110T2028_Almach.bmp';
output_plots = false;

%function [highOrderSurface] = wavefrontReconstruction(im_filename,output_plots, zernikeTerms)

% Read the image.
imageColor = imread(file_name);
image = rgb2gray(imageColor);

if output_plots
    figure, imshow(image,[]);
end

% Get parameters of the regular grid according to known bright spots.
[referenceX, referenceY, magnification] = GetGrid(image, output_plots);
% Calculate the quiver.
[arrows, idealCoords, ~] = GetQuiver(image, referenceX, referenceY, magnification, output_plots);

% Get the regular slope maps.
slopeMagnification = 1; lateralMagnification = 1;
integrationStep = 1;
[regularSlopeX,regularSlopeY,xCoordinates,yCoordinates] = ...
    Quiver2RegularSlope(arrows,idealCoords,slopeMagnification, ...
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
[removedSurface,fitCoeff] = RemoveLowOrderZernike(xCentered, yCentered, shapeDiff, zernikeTerms, output_plots);
