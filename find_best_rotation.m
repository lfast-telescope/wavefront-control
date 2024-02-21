function [best_angle] = find_best_rotation(imCrop,angle_start,angle_range)
    %This is such a stupid way to find the best rotation angle of SH but
    %it's better than wasting more time learning how to do optimization in
    %MATLAB. #afternoonignominy

    %History - 2/19/24 Warren Foster
    
    num_trials = 3;
    pk_holder = zeros(1,num_trials);
    a = linspace(angle_start-angle_range,angle_start+angle_range,num_trials);
    i=1;
    for angle = a
        rotImage = imrotate(imCrop,angle);
        x = mean(rotImage,1);
        x_smooth = medfilt1(x,5);
        [pks, locs] = findpeaks(x_smooth,'MinPeakDistance',20);
        peak_list = sort(pks,'descend');
        peak_avg = mean(peak_list(1:20));
        pk_holder(i) = peak_avg;
        i=i+1;
    end
    index = find(pk_holder == max(pk_holder));
    best_angle = a(index);
end