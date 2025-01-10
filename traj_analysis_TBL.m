clear;clc;close all

% Set path were functions will be read from
addpath(genpath('/Users/fcb/Documents/GitHub/Particles-Over-Canopy-WT/'));

fname = 'Particle_Het_Tooth_Tooth';

folderin = '/Users/fcb/Library/CloudStorage/GoogleDrive-facundo@pdx.edu/My Drive/Canopy Experiment/outputs_PTV/Toothpick/';
folderout = ['/Users/fcb/Library/CloudStorage/GoogleDrive-facundo@pdx.edu/My Drive/Canopy Experiment/analyses/' ...
    'Toothpick/' fname '/turbulent_boundary_layer/'];
mkdir(folderout)
cd(folderout)

Fs=4e3; % Frame rate

%%% Values to Normalize
L = 30; % canopy height: 30 mm
V = 1.5e3; % free-stream vel: 1500 mm/s

%% Concatenate data

if pi==pi
    trajs_conc = [];
    
    load([folderin filesep 'trajsf_' fname '_C01.mat'])
    trajs_conc = [trajs_conc tracklong];
    
    load([folderin filesep 'trajsf_' fname '_C02.mat'])
    trajs_conc = [trajs_conc tracklong];

    load([folderin filesep 'trajsf_' fname '_C03.mat'])
    trajs_conc = [trajs_conc tracklong];
    

clear tracklong 
Ine=find(arrayfun(@(X)(~isempty(X.Vx)),trajs_conc)==1);
trajs_conc = trajs_conc(Ine);
save([folderin filesep 'traj_conc_' fname],'trajs_conc','Ine','-v7.3')
end

load([folderin filesep 'traj_conc_' fname])

Ine=find(arrayfun(@(X)(~isempty(X.Vx)),trajs_conc)==1);
trajs_conc = trajs_conc(Ine);

%%
mycolormap = mycolor('#063970','#e28743');%('#063970','#eeeee4','#e28743')
color3 = [mycolormap(1,:);mycolormap((size(mycolormap,1)+1)/2,:);mycolormap(end,:)];
color1 = '#476d76';

%% Get rid of background manually 
disp('Do you want to do this part?')
%%%
% 
% for i=1:numel(trajs_conc)
% plot3(trajs_conc(i).Xf,trajs_conc(i).Yf,trajs_conc(i).Zf,'.- ');
% end
% 
% xlabel('X (streamwise)')
% ylabel('Y (vertical - antigravity)')
% zlabel('Z (spanwise)')
%%%

trajs_conc_nob = [];

for i=1:numel(trajs_conc)

    if trajs_conc(i).Xf>-17
            i/numel(trajs_conc)
           trajs_conc_nob=vertcat(trajs_conc_nob,trajs_conc(i));
    end
end

trajs_conc = trajs_conc_nob;

%% eliminate elements with no velocity
trajs_conc_nonan = [];

for i=1:numel(trajs_conc)
    if ~isnan(trajs_conc(i).Vx)
         i/numel(trajs_conc)
           trajs_conc_nonan=vertcat(trajs_conc_nonan,trajs_conc(i));
    end
end

trajs_conc = trajs_conc_nonan;

%% Rotate data to have the vertical on Y and the streamwise direction in X

for i=1:numel(trajs_conc)
    % positions
    trajs_conc_new_axis(i).Xf = trajs_conc(i).Yf; % streamwise
    trajs_conc_new_axis(i).Yf = sqrt(2)/2*(trajs_conc(i).Xf + trajs_conc(i).Zf); % vertical
    trajs_conc_new_axis(i).Zf = sqrt(2)/2*(-trajs_conc(i).Xf + trajs_conc(i).Zf); % spanwise
    % velocities
    trajs_conc_new_axis(i).Vx = trajs_conc(i).Vy; % streamwise
    trajs_conc_new_axis(i).Vy = sqrt(2)/2*(trajs_conc(i).Vx + trajs_conc(i).Vz); % vertical
    trajs_conc_new_axis(i).Vz = sqrt(2)/2*(-trajs_conc(i).Vx + trajs_conc(i).Vz); % spanwise
end

%% Plot 3D Trajectories

figure(100); clf; hold on; grid on; box on

%streamw_vel_mean = mean(vertcat(trajs_conc_new_axis.Vx));
streamw_vel_mean = mean(vertcat(trajs_conc.Vx));

for i=1:numel(trajs_conc)
i
    scatter3(trajs_conc_new_axis(i).Xf,trajs_conc_new_axis(i).Yf,trajs_conc_new_axis(i).Zf, 10, (trajs_conc_new_axis(i).Vx - streamw_vel_mean)./V, 'filled');
    %scatter3(trajs_conc(i).Xf,trajs_conc(i).Yf,trajs_conc(i).Zf, 10, (trajs_conc(i).Vx - streamw_vel_mean)./V, 'filled');
    
end
xlabel('X/L')
ylabel('Y/L')
zlabel('Z normalized (spanwise)')
colorbar

cb = colorbar(); 
ylabel(cb,'(Vx - U)/U','FontSize',20,'Rotation',270)
 
xlim([-30 30])

%stop

savefig_FC('traj3D',8,6,'pdf')
savefig_FC('traj3D',8,6,'fig')
%% Plot Position Histograms

figure(2); clf; hold on; grid on; box on

    % histogram(vertcat(trajs_conc.Xf),'FaceColor','r')
    % histogram(vertcat(trajs_conc.Yf),'FaceColor','g')
    % histogram(vertcat(trajs_conc.Zf),'FaceColor','b')

    histogram(vertcat(trajs_conc_new_axis.Xf)./L,'FaceColor',color3(1,:))
    histogram(vertcat(trajs_conc_new_axis.Yf)./L,'FaceColor',color3(2,:))
    histogram(vertcat(trajs_conc_new_axis.Zf)./L,'FaceColor',color3(3,:))

    legend({'X (normalized)','Y','Z', 'Streamw', 'Vertical'})

savefig_FC('histogram_pos',8,6,'pdf')
savefig_FC('histogram_pos',8,6,'fig')
%% Plot Velocity Histograms

figure(3); clf; hold on; grid on; box on

    % histogram(vertcat(trajs_conc.Vx),'FaceColor','r')
    % histogram(vertcat(trajs_conc.Vy),'FaceColor','g')
    % histogram(vertcat(trajs_conc.Vz),'FaceColor','b')

    histogram(vertcat(trajs_conc_new_axis.Vx)./V,'FaceColor',color3(1,:))
    histogram(vertcat(trajs_conc_new_axis.Vy)./V,'FaceColor',color3(2,:))
    histogram(vertcat(trajs_conc_new_axis.Vz)./V,'FaceColor',color3(3,:))

    legend({'Vx (normalized)','Vy','Vz', 'Streamw', 'Vertical'})

savefig_FC('histogram_vel',8,6,'pdf')
savefig_FC('histogram_vel',8,6,'fig')
%% Vertical velocity versus height bin 

figure(4); clf; hold on; grid on; box on

%streamw_vel_mean = mean(vertcat(trajs_conc_new_axis.Vx));
streamw_vel_mean = 0;

ybins = (-15:5:20);

vels_in_each_bin = cell(1, length(ybins)-1);
meanvel_in_each_bin = zeros(1, length(ybins)-1);


meanVy = cell(1, length(ybins)-1);
countVy = cell(1, length(ybins)-1);

VyBins = cell(1, length(ybins)-1);

for i = 1:numel(trajs_conc_new_axis)
    %%%%% Quantity to plot
    Yf = trajs_conc_new_axis(i).Yf;
    Vy = trajs_conc_new_axis(i).Vy;
    %%%%%

    [~, binIdx] = histc(Yf, ybins);
    
    binIdx(binIdx == 0) = 1;
    binIdx(binIdx == length(ybins)) = length(ybins) - 1;
    
    for j = 1:length(ybins)-1
        %%%% Decide if removing mean vel or not.
        VyBins{j} = [VyBins{j}; Vy(binIdx == j) - streamw_vel_mean];
        %%%%
    end
end

for j = 1:length(ybins)-1
    meanVy{j} = mean(VyBins{j});
    countVy{j} = numel(VyBins{j});
end

meanVy = cell2mat(meanVy);
countVy = cell2mat(countVy);

bar(ybins(1:end-1)./L, meanVy./V);
hold on;
errorbar(ybins(1:end-1)./L, meanVy./V, std(VyBins{j})./sqrt(countVy)./V, 'r.', 'LineWidth', 2);
xlabel('y (norm)');
ylabel('Mean Vert Vel');
title('Mean Value of Vert Vel (normalized)');

savefig_FC('vertical_vel_versus_height',8,6,'pdf')
savefig_FC('vertical_vel_versus_height',8,6,'fig')


%% Streamwise velocity versus height bin 

figure(5); clf; hold on; grid on; box on

streamw_vel_mean = mean(vertcat(trajs_conc_new_axis.Vx));

ybins = (-15:5:20);

vels_in_each_bin = cell(1, length(ybins)-1);
meanvel_in_each_bin = zeros(1, length(ybins)-1);


meanVy = cell(1, length(ybins)-1);
countVy = cell(1, length(ybins)-1);

VyBins = cell(1, length(ybins)-1);

for i = 1:numel(trajs_conc_new_axis)
    %%%%% Quantity to plot
    Yf = trajs_conc_new_axis(i).Yf;
    Vy = trajs_conc_new_axis(i).Vx;
    %%%%%

    [~, binIdx] = histc(Yf, ybins);
    
    binIdx(binIdx == 0) = 1;
    binIdx(binIdx == length(ybins)) = length(ybins) - 1;
    
    for j = 1:length(ybins)-1
        %%%% Decide if removing mean vel or not.
        VyBins{j} = [VyBins{j}; Vy(binIdx == j) - streamw_vel_mean];
        %%%%
    end
end

for j = 1:length(ybins)-1
    meanVy{j} = mean(VyBins{j});
    countVy{j} = numel(VyBins{j});
end

meanVy = cell2mat(meanVy);
countVy = cell2mat(countVy);

bar(ybins(1:end-1)./L, meanVy./V);
hold on;
errorbar(ybins(1:end-1)./L, meanVy./V, std(VyBins{j})./sqrt(countVy)./V, 'r.', 'LineWidth', 2);
xlabel('y (norm)');
ylabel('Mean Streamw Vel ');
title('Mean Value of Streamw Vel (Norm)');

savefig_FC('streamwise_vel_versus_height',8,6,'pdf')
savefig_FC('streamwise_vel_versus_height',8,6,'fig')

stop
%% Gradient (over height) of vertical velocity 

figure(6); clf; hold on; grid on; box on


for i=1:numel(trajs_conc_new_axis)
    
    y = trajs_conc_new_axis(i).Yf;
    
    gradients(i).dVx_dy = diff(trajs_conc_new_axis(i).Vx) ./ diff(y);
    gradients(i).dVy_dy = diff(trajs_conc_new_axis(i).Vy) ./ diff(y);
    
end

plot(vertcat(gradients.dVx_dy),'.')


















