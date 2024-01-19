function [regularSlopeX,regularSlopeY,xCoordinates,yCoordinates] = Quiver2RegularSlope(arrows,idealCoords,slopeMagnification,lateralMagnification,integrationStep)

% This function is used to get the regular slope maps based on irregular
% slope maps.
% INPUTS:
%   arrows: N x 2 matrix showing displacements on the image;
%   idealCoords: N x 2 matrix showing x and y coordinates of available
% ideal points;
%   slopeMagnification: magnification used to turn arrows to slope data;
%   lateralMagnification: magnification of the mirror, which indicates the 
% actual spacing on the mirror coverd by a single camera pixel (mm);
%   integrationStep: integration step, unit in mm.
% OUTPUTS: 
%   regularSlopeX/regularSlopeY: regular x/y slope map;
% HISTORY:
%   2023-11-03 - Yiyang Huang - rewrite from 'RegularizeSlope'

% Decompose data.
slope = slopeMagnification * arrows;
mirrorPos = lateralMagnification * idealCoords;

% Get regular grids showing x and y coordinates of to-be-fitted positions.
xStart = round(min(mirrorPos(:,1))); xEnd = round(max(mirrorPos(:,1)));
yStart = round(min(mirrorPos(:,2))); yEnd = round(max(mirrorPos(:,2)));
rangeX = xStart : integrationStep : xEnd;
rangeY = yStart : integrationStep : yEnd;
[xCoordinates,yCoordinates] = meshgrid(rangeX,rangeY); % give x and y coordinates of fitted points

% Fit slope maps.
regularSlopeX = griddata(mirrorPos(:,1),mirrorPos(:,2),slope(:,1), ...
    xCoordinates,yCoordinates);
regularSlopeY = griddata(mirrorPos(:,1),mirrorPos(:,2),slope(:,2), ...
    xCoordinates,yCoordinates);
regularSlopeX = integrationStep * regularSlopeX; % unify integration spacing
regularSlopeY = integrationStep * regularSlopeY;

% Display and save regular slope maps.
% figure, mesh(xCoordinates,yCoordinates,regularSlopeX); 
% view([0 90]); axis equal
% figure, mesh(xCoordinates,yCoordinates,regularSlopeY); 
% view([0 90]); axis equal
% figure, imhandle = meshc(xCoordinates,yCoordinates,regularSlopeX); 
% view([0 90]); axis equal
% imhandle(2).EdgeColor = 'k'; imhandle(2).ZLocation = 'zmax';
% % imhandle(2).LevelStep = 5e-4; 
% caxis([-2e-3 4.5e-3]); % (!)
% colormap('jet'); colorbar; title('Regular slope in x direction');
% exportgraphics(gca,['regular slope X','.png']);
% figure, imhandle = meshc(xCoordinates,yCoordinates,regularSlopeY); 
% view([0 90]); axis equal
% imhandle(2).EdgeColor = 'k'; imhandle(2).ZLocation = 'zmax';
% % imhandle(2).LevelStep = 5e-4; 
% caxis([-1e-3 3.5e-3]); % (!)
% colormap('jet'); colorbar; title('Regular slope in y direction');
% exportgraphics(gca,['regular slope Y','.png']);

end
