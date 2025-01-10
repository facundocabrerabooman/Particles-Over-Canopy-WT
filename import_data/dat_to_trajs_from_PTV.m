clear;clc;close all

% Set path were functions will be read from
addpath(genpath('/Users/fcb/Documents/GitHub/Particles-Over-Canopy-WT/'));

fname = 'Particle_Het_Tooth_Tooth_C02';

folderin = '/Users/fcb/Library/CloudStorage/GoogleDrive-facundo@pdx.edu/My Drive/Canopy Experiment/outputs_PTV/Toothpick';
folderout = folderin;
cd(folderin)

Fs=4e3; % Frame rate

%% Import data

d = dat_to_mat(folderin, fname);
%save(fname,'d','-v7.3')

%% Convert to structure

tracklong.Xf = d(:,2);
tracklong.Yf = d(:,3);
tracklong.Zf = d(:,4);

tracklong.Vx = d(:,5);
tracklong.Vy = d(:,6);
tracklong.Vz = d(:,7);

tracklong.Ax = d(:,8);
tracklong.Ay = d(:,9);
tracklong.Az = d(:,10);

stop
save(['trajsf_' fname '.mat'],'Ine','tracklong','-v7.3')


