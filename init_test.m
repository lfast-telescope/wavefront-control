path = '/home/warrenbfoster/Documents/on-sky/20241204/';
list_of_subfolders = dir(path);
chosen_subfolder = list_of_subfolders(4).name;

subfolder_path = append(path,chosen_subfolder,'/');
list_of_files = {dir(subfolder_path).name};

zernikeTerms = [1,2,3,4];

im_filename = append(subfolder_path,list_of_files(4));
output_plots = false;
filename = append(path, 'reconstruction.mat');

[highOrderSurface, xCoordinates, yCoordinates] = wavefrontReconstruction(im_filename,output_plots,zernikeTerms, filename);


%x,y,z = wavefrontReconstruction()