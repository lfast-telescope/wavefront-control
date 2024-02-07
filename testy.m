clearvars; clc; close all;

%Testbench for "wavefrontReconstruction.m"

output_plots = false;
zernikeTerms = [1,2,3,4];

folder_path = 'C:\Users\warrenbfoster\Documents\LFAST\on-sky\20240118\';
file_list = dir(folder_path);

for iteration = 3:10
    im_filename = file_list(iteration).name;
    [highOrderSurface,xCoordinates,yCoordinates] = wavefrontReconstruction(append(folder_path,im_filename),output_plots, zernikeTerms);
    
    sortedVals = sort(highOrderSurface);
    lowThresh = sortedVals(uint32(numel(highOrderSurface)*0.01));
    highThresh = sortedVals(uint32(numel(highOrderSurface)*0.99));
    figure, imhandle = meshc(xCoordinates,yCoordinates,highOrderSurface);
    
    imhandle(2).EdgeColor = 'k'; imhandle(2).ZLocation = 'zmax';
    imhandle(2).LevelStep = 100; view([0 90]); 
    axis equal; title(im_filename); colormap('jet'); colorbar;
    hold off

end