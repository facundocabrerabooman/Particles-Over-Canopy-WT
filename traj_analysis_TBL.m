clear;clc;close all

% Set path were functions will be read from
addpath(genpath('/Users/fcb/Documents/GitHub/Particles-Over-Canopy-WT/'));

fname = 'Particle_Het_Foam_Gap';

folderin = '/Users/fcb/Library/CloudStorage/GoogleDrive-facundo@pdx.edu/My Drive/Canopy Experiment/outputs_PTV/Foam/';
folderout = ['/Users/fcb/Library/CloudStorage/GoogleDrive-facundo@pdx.edu/My Drive/Canopy Experiment/analyses/' ...
    'Foam/Particle_Het_Foam_Gap/turbulent_boundary_layer/'];
mkdir(folderout)
cd(folderout)

Fs=4e3; % Frame rate

%% Concatenate data

if pi==pi
    trajs_conc = [];
    
    load([folderin filesep 'trajsf_Particle_Het_Foam_Gap_C01.mat'])
    trajs_conc = [trajs_conc tracklong];
    
    load([folderin filesep 'trajsf_Particle_Het_Foam_Gap_C02.mat'])
    trajs_conc = [trajs_conc tracklong];
    
    load([folderin filesep 'trajsf_Particle_Het_Foam_Gap_C03.mat'])
    trajs_conc = [trajs_conc tracklong];
    
clear tracklong 
Ine=find(arrayfun(@(X)(~isempty(X.Vx)),trajs_conc)==1);
trajs_conc = trajs_conc(Ine);
save([folderin filesep 'traj_conc_Particle_Het_Foam_Gap'],'trajs_conc','Ine','-v7.3')
end

load([folderin filesep 'traj_conc_' fname])

Ine=find(arrayfun(@(X)(~isempty(X.Vx)),trajs_conc)==1);
trajs_conc = trajs_conc(Ine);
%%
mycolormap = mycolor('#063970','#e28743');%('#063970','#eeeee4','#e28743')
color3 = [mycolormap(1,:);mycolormap((size(mycolormap,1)+1)/2,:);mycolormap(end,:)];
color1 = '#476d76';
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

figure(10); clf; hold on; grid on; box on

%streamw_vel_mean = mean(vertcat(trajs_conc_new_axis.Vx));
streamw_vel_mean = mean(vertcat(trajs_conc.Vx));

for i=1:numel(trajs_conc_new_axis)

    %scatter3(trajs_conc_new_axis(i).Xf,trajs_conc_new_axis(i).Yf,trajs_conc_new_axis(i).Zf, 10, trajs_conc_new_axis(i).Vx - streamw_vel_mean, 'filled');
    scatter3(trajs_conc(i).Xf,trajs_conc(i).Yf,trajs_conc(i).Zf, 10, trajs_conc(i).Vx - streamw_vel_mean, 'filled');

end
xlabel('X (streamwise)')
ylabel('Y (vertical - antigravity)')
zlabel('Z (spanwise)')
 
colorbar

savefig_FC('traj3D',8,6,'pdf')
savefig_FC('traj3D',8,6,'fig')
%% Plot Position Histograms

figure(2); clf; hold on; grid on; box on

    % histogram(vertcat(trajs_conc.Xf),'FaceColor','r')
    % histogram(vertcat(trajs_conc.Yf),'FaceColor','g')
    % histogram(vertcat(trajs_conc.Zf),'FaceColor','b')

    histogram(vertcat(trajs_conc_new_axis.Xf),'FaceColor',color3(1,:))
    histogram(vertcat(trajs_conc_new_axis.Yf),'FaceColor',color3(2,:))
    histogram(vertcat(trajs_conc_new_axis.Zf),'FaceColor',color3(3,:))

    legend({'x','y','z', 'Streamw', 'Vertical'})

savefig_FC('histogram_pos',8,6,'pdf')
savefig_FC('histogram_pos',8,6,'fig')
%% Plot Velocity Histograms

figure(3); clf; hold on; grid on; box on

    % histogram(vertcat(trajs_conc.Vx),'FaceColor','r')
    % histogram(vertcat(trajs_conc.Vy),'FaceColor','g')
    % histogram(vertcat(trajs_conc.Vz),'FaceColor','b')

    histogram(vertcat(trajs_conc_new_axis.Vx),'FaceColor',color3(1,:))
    histogram(vertcat(trajs_conc_new_axis.Vy),'FaceColor',color3(2,:))
    histogram(vertcat(trajs_conc_new_axis.Vz),'FaceColor',color3(3,:))

    legend({'x','y','z', 'Streamw', 'Vertical'})

savefig_FC('histogram_vel',8,6,'pdf')
savefig_FC('histogram_vel',8,6,'fig')
%% Vertical velocity versus height bin 

figure(4); clf; hold on; grid on; box on

%streamw_vel_mean = mean(vertcat(trajs_conc_new_axis.Vx));
streamw_vel_mean = 0;

ybins = (-25:5:20);

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

bar(ybins(1:end-1), meanVy);
hold on;
errorbar(ybins(1:end-1), meanVy, std(VyBins{j})./sqrt(countVy), 'r.', 'LineWidth', 2);
xlabel('y (mm)');
ylabel('Mean Vert Vel');
title('Mean Value of Vert Vel');

savefig_FC('vertical_vel_versus_height',8,6,'pdf')
savefig_FC('vertical_vel_versus_height',8,6,'fig')


%% Streamwise velocity versus height bin 

figure(5); clf; hold on; grid on; box on

streamw_vel_mean = mean(vertcat(trajs_conc_new_axis.Vx));

ybins = (-25:5:20);

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

bar(ybins(1:end-1), meanVy);
hold on;
errorbar(ybins(1:end-1), meanVy, std(VyBins{j})./sqrt(countVy), 'r.', 'LineWidth', 2);
xlabel('y (mm)');
ylabel('Mean Streamw Vel');
title('Mean Value of Streamw Vel');

savefig_FC('streamwise_vel_versus_height',8,6,'pdf')
savefig_FC('streamwise_vel_versus_height',8,6,'fig')

%% Gradient (over height) of vertical velocity 

figure(6); clf; hold on; grid on; box on


for i=1:numel(trajs_conc_new_axis)
    
    y = trajs_conc_new_axis(i).Yf;
    
    gradients(i).dVx_dy = diff(trajs_conc_new_axis(i).Vx) ./ diff(y);
    gradients(i).dVy_dy = diff(trajs_conc_new_axis(i).Vy) ./ diff(y);
    
end

plot(vertcat(gradients.dVx_dy),'.')


















