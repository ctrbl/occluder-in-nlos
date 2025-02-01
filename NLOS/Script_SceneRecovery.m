%% Test script recovery from mitsuba rendered images
% Team 2: Note that your task is to explore the effect of the hidden
% occluder's position on the reconstruction accuracy. To achieve this:
% 1) Add noise to the simulated measurement (e.g., using awgn(y_meas_r,SNR,'measured') function)
%    for a fixed SNR, say SNR = 30.
% 2) Recover the hidden scene using the noisy measurement
% 3) Change occluder location (in 'load_expt_config_measurements.m' script and 
%    repeat from step 1.

clear variables;
close all

clc;
title('Scene')

% DISCRETIZATION
Ndiscr_mon = 6; % Higher means finer scene block discretization (6 is usually ok!)

load_expt_config_measurements

downsamp_factor = 4;
%(numPixels/2^downsamp_factor num of pixels)
% downsamp_factor = 2 usually a good choice

tic
[ simA ] = SimForwardModelOnly(simuParams,downsamp_factor);
toc


%%

figure(1); clf;

% Make a synthetic hidden scene and display it along with its measurements
% I give a simple scene here, but modify it to get a more interesting scene
% (e.g., you could make a letter or number or plus sign)
% Red channel
sample_scene_r = zeros(NumBlocks_sim);
sample_scene_r(1,13) = 1; % Top left pixel
sample_scene_r(1,1) = 1; % Top right pixel
% Green channel
sample_scene_g = zeros(NumBlocks_sim);
sample_scene_g(7,13) = 1; % Bottom left pixel
sample_scene_g(1,1) = 1; % Top right pixel
% Blue channel
sample_scene_b = zeros(NumBlocks_sim);
sample_scene_b(7,13) = 0; % Bottom left pixel
sample_scene_b(1,1) = 1; % Top right pixel


% This is just for displaying a plot of the scene
MonitorPattern_comp_r = fliplr(kron(sample_scene_r,ones(Ndiscr_mon,Ndiscr_mon)));
MonitorPattern_comp_g = fliplr(kron(sample_scene_g,ones(Ndiscr_mon,Ndiscr_mon)));
MonitorPattern_comp_b = fliplr(kron(sample_scene_b,ones(Ndiscr_mon,Ndiscr_mon)));

MonitorPattern_comp = cat(3,MonitorPattern_comp_r, MonitorPattern_comp_g, MonitorPattern_comp_b);

subplot(121)
imagesc(MonitorPattern_comp,[0 1]);  axis equal tight
title('Scene')
% Remember equal amounts of Red + Green gives yellow.
% Remember equal amounts of Red + Green + Blue gives white.




% Simulate noiseless measurement
numSimPixMeas = numPixels/(2^downsamp_factor);
y_meas_display = cat(3, reshape((simA*sample_scene_r(:)), fliplr(numSimPixMeas)),...
                reshape((simA*sample_scene_g(:)), fliplr(numSimPixMeas)),...
                reshape((simA*sample_scene_b(:)), fliplr(numSimPixMeas)));
y_meas_simulated = cat(3, reshape((simA*sample_scene_r(:)), (numSimPixMeas)),...
                reshape((simA*sample_scene_g(:)), (numSimPixMeas)),...
                reshape((simA*sample_scene_b(:)), (numSimPixMeas)));

subplot(122);
imagesc(y_meas_display./max(y_meas_display(:)))
axis square
colormap(gca, 'gray')
title('Measurement')


%% Pseudo-Inverse inverse (Real data)
% This reconstruction uses a real photograph of a setup in my Lab (ENG
% 204) you may set up a meeting to stop by if you want to see the setup
% and/or collect your own dataset.

renderTestIm = 'RGB_bars_crop';
[y_meas_real] = LoadRenderPhoto(renderTestIm, downsamp_factor,numPixels);

sr = 1;
sg = 1;
sb = 1;



figure(3)

y_meas_r = y_meas_real(:,:,1);
y_meas_g = y_meas_real(:,:,2);
y_meas_b = y_meas_real(:,:,3);
scene_est = cat(3, reshape(simA\y_meas_r(:),NumBlocks_sim),...
                   reshape(simA\y_meas_g(:),NumBlocks_sim),...
                   reshape(simA\y_meas_b(:),NumBlocks_sim));
% This medfilt2 below is used to supress artifacts in reconstruction (not
% needed for simulated noiseless data)
final_im_m1a = cat(3,medfilt2(scene_est(:,:,1)),medfilt2(scene_est(:,:,2)),medfilt2(scene_est(:,:,3)));
imagesc(fliplr(final_im_m1a)*180)
title('Reconstruction - Real Data')
axis equal tight



%% Pseudo-Inverse Simulated Data


figure(4)

y_meas_r = y_meas_simulated(:,:,1);
y_meas_g = y_meas_simulated(:,:,2);
y_meas_b = y_meas_simulated(:,:,3);

% This is where you want to add noise to each of the three color channels

scene_est = cat(3, reshape(simA\y_meas_r(:),NumBlocks_sim),...
                   reshape(simA\y_meas_g(:),NumBlocks_sim),...
                   reshape(simA\y_meas_b(:),NumBlocks_sim));
imagesc(scene_est)
title('Reconstruction - Simulated Data')
axis equal tight

