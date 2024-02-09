%% Split tiff images into one per camera
close all
clear
clc

addpath(genpath('/Users/fcb/Documents/GitHub/Particles-Over-Canopy-WT'));

% input
fpath0 = '/Volumes/Canopy_1/10PPI_C01/';
fpath = '/Volumes/Canopy_1/10PPI_C01_split/';
fpathout = '/Volumes/Canopy_1/grayscale/';
mkdir(fpath); mkdir(fpathout);
mkdir([fpath filesep 'Camera1']);mkdir([fpath filesep 'Camera2']); % create folders for spit images
mkdir([fpathout filesep 'Camera1']);mkdir([fpathout filesep 'Camera2']); % create folders for postprocessed spit images


image_list = dir([fpath0 filesep '*.tif']);
img_num = size(image_list,1);


for i=1:img_num

    i/img_num
    fname=[fpath0 filesep 'B' num2str(i,'%05d') '.tif'];
    Im_original = imread(fname);
    Imc1 = Im_original(1:960,:);
    Imc2 = Im_original(961:end,:);

    imwrite(uint16(Imc1),[fpath filesep 'Camera1' filesep 'B' num2str(i,'%05d') '.tif'])
    imwrite(uint16(Imc2),[fpath filesep 'Camera2' filesep 'B' num2str(i,'%05d') '.tif'])

end