close all
clear
clc

part_radius = 3;
addpath(genpath('/Users/fcb/Documents/GitHub/Particles-Over-Canopy-WT'));

% input
fpath0 = '/Volumes/Canopy_1/10PPI_C01/';
fpath = '/Volumes/Canopy_1/10PPI_C01_split/';
fpathout = '/Volumes/Canopy_1/grayscale/';
mkdir(fpath); mkdir(fpathout);
mkdir([fpath filesep 'Camera1']);mkdir([fpath filesep 'Camera2']); % create folders for spit images
mkdir([fpathout filesep 'Camera1']);mkdir([fpathout filesep 'Camera2']); % create folders for postprocessed spit images

image_list = dir([fpath filesep 'Camera1' filesep '*.tif']);
Iend = size(image_list,1); 
Istart=1;
Iend=200;

%%% get Backgrounds
%cam 1
bkg1_originalSize = getBkg(fpath,'Camera1',Istart,Iend,100,[]);
%cam 2
bkg2_originalSize = getBkg(fpath,'Camera2',Istart,Iend,100,[]);

counter = 0;
for k=Istart:Iend
    counter = counter+1;
    k/Iend

    %%% cam1
    fname=[fpath filesep 'Camera1' filesep 'B' num2str(k,'%05d') '.tif'];
    Im1_originalSize = imread(fname);

    intensity_thr = 40;
    Imc1 = Ryan_preprocessing(Im1_originalSize,intensity_thr,part_radius,1,bkg1_originalSize);

    fnameo = ['cam1_' num2str(counter,'%05d') '.tif'];
    imwrite(uint16(Imc1),[fpathout filesep 'Camera1' filesep fnameo])
    %                 figure(10);
    %                 subplot(2,1,1);imagesc(Im1_originalSize);axis equal
    %                 subplot(2,1,2);imagesc(Im1p);axis equal
    %                 pause(0.05)


    %%% cam2
    fname=[fpath filesep 'Camera2' filesep 'B' num2str(k,'%05d') '.tif'];
    Im2_originalSize = imread(fname);

    intensity_thr = 25;
    Imc2 = Ryan_preprocessing(Im2_originalSize,intensity_thr,part_radius,1.5,bkg2_originalSize);

    fnameo = ['cam2_' num2str(counter,'%05d') '.tif'];
    imwrite(uint16(Imc2),[fpathout filesep 'Camera2' filesep fnameo]);

    %                figure(10);
    %                  subplot(2,1,1);imagesc(Im2_originalSize);axis equal
    %                  subplot(2,1,2);imagesc(Im2p);axis equal
    %                  pause(0.05)



end
stop
%% Test camera 1

if 1==pi
    figure;
    subplot(2,1,1)
    imshow(Im1_originalSize);axis equal;
    [C, ~] = imfindcircles(imbinarize(Im1p),[1 15]);
    viscircles(C,ones(size(C,1),1).*12,'Color','r')


    subplot(2,1,2)
    imagesc(Im1p);axis equal; hold on
    viscircles(C,ones(size(C,1),1).*12,'Color','r')

    figure;imagesc(Im1t);axis equal;hold on
    %viscircles(C,ones(size(C,1),1).*5,'Color','r')

end
%% Test camera 2
if 1==pi
    figure;
    subplot(2,1,1)
    imshow(Im2_originalSize);axis equal;
    [C, ~] = imfindcircles(imbinarize(Im2p),[1 15]);
    viscircles(C,ones(size(C,1),1).*12,'Color','r');


    subplot(2,1,2)
    imagesc(Im2p);axis equal; hold on
    viscircles(C,ones(size(C,1),1).*12,'Color','r');

    figure;imagesc(Im2t);axis equal;hold on

end
%% Test camera 3
if 1==pi
    figure;
    subplot(2,1,1)
    imshow(Im3_originalSize);axis equal;
    [C, ~] = imfindcircles(imbinarize(Im3p),[1 15]);
    viscircles(C,ones(size(C,1),1).*12,'Color','r');


    subplot(2,1,2)
    imagesc(Im3p);axis equal; hold on
    viscircles(C,ones(size(C,1),1).*12,'Color','r');

    figure;imagesc(Im3t);axis equal;hold on
end