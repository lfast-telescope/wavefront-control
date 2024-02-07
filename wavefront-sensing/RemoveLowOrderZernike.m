function [highOrderMap,fitCoeff] = RemoveLowOrderZernike(xCoordinates, yCoordinates, surface, zernikeTerms, output_plots)


% This function removes some low order Zernike terms.
% INPUT:
%   xCoordinates: x coordinates of the surface;
%   yCoordinates: y coordinates of the surface;
%   surface: surface shape of the measured mirror;
%   zernikeTerms: the Zernike terms to be removed.
% OUTPUT:
%   removedSurface: the surface after removing the low order surface.
% HISTORY:
%   2023-07-25 - Yiyang Huang - initial implementation
%   2023-11-09 - Yiyang Huang - adapt to LFAST
%   2024-01-12 - Warren Foster - enable plot suppression
%   2024-01-15 - Yiyang Huang - encircle the low order map

% Default Zernike terms to be removed.
if nargin < 4
    zernikeTerms = [1,2,3,4,5,6];
end


% Normalize x and y coordinates into a unit circle.
[theta, rho] = cart2pol(xCoordinates,yCoordinates);
abRho = abs(rho); radius = max(abRho(:));
rho = rho / radius;
zernikeVal = Zernike(rho, theta, zernikeTerms); % Zernike values with different orders


% Get coefficients of low order Zernike terms. 
surface0 = surface(:); 
termsNum = length(zernikeTerms);
zernikeVal0 = [];
for i = 1:termsNum
    currentZernike = zernikeVal(:,:,i);
     % figure(100), imagesc(currentZernike); 
     % colormap('jet'); colorbar; axis equal
     % pause(2); close(100);
    zernikeVal0 = [zernikeVal0 currentZernike(:)];
end
% Remove nan data since they'll mess up fitting.
badData = any(isnan(zernikeVal0),2) | any(isnan(surface0),2);
zernikeVal0 = zernikeVal0(~badData,:);
surface0 = surface0(~badData,:);
% Solve equation of form: surface0 = zernikeVal0*fitCoeff 
% (lowOrderSurface = ZernikeTerms*ZernikeCoefficients)
fitCoeff = zernikeVal0\surface0;


% Fit the low-order surface and subtract it from the original image.
lowOrderMap = zeros(size(surface));
lowOrderInfo = 'Zernike Coefficients: ';
for i = 1:termsNum
    lowOrderMap = lowOrderMap + fitCoeff(i) * zernikeVal(:,:,i);
    lowOrderInfo = [lowOrderInfo,num2str(i),' ~ ',num2str(fitCoeff(i)),'; '];
end

highOrderMap = surface - lowOrderMap;

if output_plots
    figure, mesh(xCoordinates,yCoordinates,lowOrderMap); view([0 90]); 
    axis equal; title('Low order map'); colormap('jet'); colorbar; 
    subtitle(lowOrderInfo,'FontSize',8);
    exportgraphics(gca,'low order surface.png');

    sortedVals = sort(surface0);
    lowThresh = sortedVals(uint32(numel(surface0)*0.01));
    highThresh = sortedVals(uint32(numel(surface0)*0.99));
    figure, imhandle = meshc(xCoordinates,yCoordinates,highOrderMap);
    
    imhandle(2).EdgeColor = 'k'; imhandle(2).ZLocation = 'zmax';
    imhandle(2).LevelStep = 100; view([0 90]); 
    axis equal; title('High order surface'); colormap('jet'); colorbar;
    %   caxis([lowThresh highThresh]);
    exportgraphics(gca,'high order surface.png');
end
end
