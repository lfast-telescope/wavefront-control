path = 'C:/Users/warrenbfoster/OneDrive - University of Arizona/Documents/LFAST/on-sky/20241204/SHWFS/';
list_of_subfolders = dir(path);

for chosen_subfolder = {list_of_subfolders.name}
    if numel(chosen_subfolder{:}) > 2
        subfolder_path = append(path,chosen_subfolder,'/');
        list_of_files = {dir(subfolder_path{:}).name};
        output_plots = false;
        
        zernikeTerms = [1,2,3,4];
        
        cell_holder = cell(size(list_of_files));
        
        i = 0;
        for file = list_of_files
            i=i+1;
            if numel(file{1})>2
                im_filename = append(subfolder_path,file); 
                if ~isfolder(im_filename)
                    filename = append(path, 'reconstruction.mat'); 
                    [highOrderSurface, xCoordinates, yCoordinates] = wavefrontReconstruction(im_filename,output_plots,zernikeTerms, filename);
                    cell_holder{i}=highOrderSurface;
                end
            end
        end
        
        mean_counter = 0;
        running_sum = zeros(size(highOrderSurface));
        for val = cell_holder
            if ~isempty(val{1})
                if all(size(val{1}) == size(running_sum))
                    running_sum = running_sum + val{1};
                    mean_counter = mean_counter + 1;
                end
            end
        end
        output_val = running_sum / mean_counter;
        
        save(append(path, 'average_', chosen_subfolder{:}, '.mat'),"output_val")
    end
end