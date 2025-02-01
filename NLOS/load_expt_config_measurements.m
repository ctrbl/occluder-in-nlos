%% This script holds the measurement configurations of the experimental dataset.
% John Murray-Bruce at USF
% Team 2: The parameters you should be changing are the entries of
% Occ_LLcorner (i.e., the occluder's position in [x y z]).

numPixels = [2272 1888]; % Number of pixels in digital camera

NumBlocks_sim = [7 13];
NumBlocks_cal = [7 13];
D = 1.249;

% OCCLUDER DEFINITION
Occluder_size = [0.1 0.1]; % [Length width] in that order

Occ_LLcorner = [.575    0.5125    0.2875]; % Position of the (lower left corner) of the square occluder in [x y z]
Occluder = [Occ_LLcorner; Occ_LLcorner + [Occluder_size(1) 0 Occluder_size(2)] ];


%CAMERA CONFIGURATION
FOV_size = [0.7 .585];
FOV_LLCorner = [1.4665-0.7 0.6/100];
FOV_cord = [FOV_LLCorner; FOV_LLCorner + FOV_size];


NumBlocks_col = NumBlocks_cal(2);
NumBlocks_row = NumBlocks_cal(1);
ScreenSize = [0.6 0.335];
IlluminationBlock_Size = (ScreenSize./[NumBlocks_col NumBlocks_row]).*[1 1];


% Set simulation of forward model parameters!
simuParams.NumBlocks = NumBlocks_sim;
simuParams.Ndiscr_mon = Ndiscr_mon;
simuParams.numPixels = numPixels;
simuParams.D = D;
simuParams.FOV_cord = FOV_cord;
simuParams.Occluder = Occluder;
simuParams.viewAngleCorrection = false;

simuParams.IlluminationBlock_Size = IlluminationBlock_Size;
simuParams.Mon_Offset = [0 0.185];
