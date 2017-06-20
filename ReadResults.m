clear;
Original_image_dir = 'Results';
fpath = fullfile(Original_image_dir, 'USPS*.mat');
mat_dir  = dir(fpath);
mat_num = length(mat_dir);


for i = 1 : mat_num
    fprintf([mat_dir(i).name '\n']);
    eval(['load Results/' mat_dir(i).name]);
end
