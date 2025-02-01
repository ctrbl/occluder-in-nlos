function [ simA, Discr ] = SimForwardModelOnly(simuParams, downsampling)
%GETLINEARFORWARDMODEL Simulates the forward operator (measurement matrix)
% for the monitor/camera experimental setup, using the set parameters
% provided in simulationParams (of type struct).
%
%   Usage:
%   [ simA, Discr ] = GetLinearForwardModel(simuParams)
%   [ simA, Discr ] = GetLinearForwardModel(simuParams,downsamp)
%   Input:
%       * simuParams (struct): containing details of experimental configuration.
%                             (1) simuParams.NumBlocks (Numv
%                             (2) simuParams.Ndiscr_mon
%                             (3) simuParams.Occluder
%                             (4) simuParams.viewAngleCorrection
%                             (5) simuParams.D
%                             (6) simuParams.FOV_cord;
%                             (7) simuParams.numPixels
%                             (8) simuParams.IlluminationBlock_Size
%       * downsampling (scalar): Douwnsampling factor by which to downsample the measurements.
%   Output:
%       * simA:           The discrete forward model matrix.
%       * Discr:          Hidden scene axes grid.

% Last Modified by $John Murray-Bruce$ at Boston University
% v1.0 05-Sept-2017 10:28:32


NumBlocks = simuParams.NumBlocks;
numPixels = simuParams.numPixels;
if length(numPixels)==1
    numPixels = [numPixels numPixels];
end
if downsampling>0
    numPixels = floor(numPixels/(2^downsampling));
end
nDiscr = simuParams.Ndiscr_mon;
D = simuParams.D;
FOV_cord = simuParams.FOV_cord;
Occluder = simuParams.Occluder;
viewAngleCorrection = simuParams.viewAngleCorrection;
IlluminationBlock_size = simuParams.IlluminationBlock_Size;
Mon_Offset = simuParams.Mon_Offset;

[rowCount, colCount] = meshgrid(1:(NumBlocks(1)),1:(NumBlocks(2)));
rowCount = rowCount'; colCount = colCount';
% Note this not so natural ordering used below, due to the file naming
% convention for the monitor blocks used by Charlie's capture code.
ActiveBlock = [rowCount(:), colCount(:)];

% Preallocate memory for array.
totalBlocks = NumBlocks(1)*NumBlocks(2);
simA = zeros(prod(numPixels),totalBlocks);


% SIMULATE CAMERA MEASUREMENTS
parfor jj=1:totalBlocks
    Im = MonitorSimulateForwardModel(NumBlocks, nDiscr, ActiveBlock(jj,:),...
        numPixels, D, FOV_cord, [Occluder],viewAngleCorrection,...
        IlluminationBlock_size, Mon_Offset);
    Im = flipud(Im);
    simA(:,jj) = Im(:);
end

Discr.Wall_xdiscr = [linspace(FOV_cord(1,1),FOV_cord(2,1),numPixels(1))];

Discr.Wall_zdiscr= linspace(FOV_cord(1,2),FOV_cord(2,2),numPixels(2));
end


function [ Image, MonPattern, Wall_discr ] = MonitorSimulateForwardModel(NumBlocks,...
    Ndiscr_mon, ActiveBlocks_coord, numPixels, Dist_MonWall, FOV, Occluder,...
    viewAngleCorrection,IlluminationBlock_size, Mon_Offset)
%MONITORSIMULATEFORWARDMODEL simulates the measurements on the visible wall
% (i.e. the camera measurements) for a specified active block (scene pixel)
% on the non-visible monitor.
%   Usage:
%   [ Image, MonPattern, Wall_discr ] = MonitorSimulateForwardModel(NumBlocks,
%    Ndiscr_mon, ActiveBlocks_coord, numPixels, Dist_MonWall, FOV, Occluder,
%    viewAngleCorrection,IlluminationBlock_size, Mon_Offset)
%   Input:
%       * NumBlocks (2-vect): Num of possible scene/monitor blocks for display.
%       * Ndiscr_mon (scalar): Num of point sources per dimension (of monitor block).
%       * ActiveBlocks_coord (2-vect): [row col] of scene/monitor pixel.
%       * numPixels (scalar): Number of camera pixels per color channel to simulate.
%       * Dist_MonWall (scalar): Distance between monitor and visible wall.
%       * FOV [m] (2-vect): Camera FoV on the visible wall ([x z]m)
%       * Occluder [m] (3-vector): Lower left corner of occluder ([x y z]).
%       * viewAngleCorrection (boolean): View angle correction On/Off.
%       * IlluminationBlock_size[m] (2-vector): Size of monitor blocks ([x,z]).
%       * Mon_Offset[m] (2-vector): Lower left corner of screen ([x, z] only).
%   Output:
%       * Image: Camera measurements (column of measurement matrix).
%       * MonPattern: Pattern on the monitor (visualizes ActiveBlocks_coord).
%       * Wall_discr: The [x z] dixretization of visible wall.

% Last Modified by $John Murray-Bruce$ at Boston University.
% v1.0: 20-Nov-2017 9:13:32 for sharing!
Occ_present = true;
if isempty(Occluder)
    Occ_present = false;
elseif Occluder==false
    Occ_present = false;
else
    occ_stickpresent = false;   % Make this true if you want the occluder stick
    occstick(:,1) = [-0.002/2 0.002/2]' + 0.5*(Occluder(1,1)+Occluder(2,1));
    occstick(:,3) = [0 Occluder(1,3)]';
    occstick(:,2) = Occluder(2,2);
end
% occ_stickpresent = false;

% DISCRETIZATION
NumBlocks_col = NumBlocks(2);
NumBlocks_row = NumBlocks(1);


% WALL DEFINITION
LambertianWall_xlim = FOV(:,1); %[0.5335 0.988];
LambertianWall_y = Dist_MonWall;
LambertianWall_zlim = FOV(:,2); %[0 0.51];
Wall_xdiscr = linspace(LambertianWall_xlim(1),LambertianWall_xlim(2),numPixels(1));
Wall_zdiscr = linspace(LambertianWall_zlim(1),LambertianWall_zlim(2),numPixels(2));


% MONITOR DEFINITION
if nargin<9
    ScreenSize = [0.337 0.271];
    IlluminationBlock_size =0.5*(ScreenSize(1)/(NumBlocks_col+0.5)+...
        (ScreenSize(2)-0.074)/NumBlocks_row);
    IlluminationBlock_size = [IlluminationBlock_size IlluminationBlock_size];
    
    MonitorOffset_x = 0.091;
    MonitorOffset_z = 0.1445;
elseif nargin==9
    MonitorOffset_x = 0.091;
    MonitorOffset_z = 0.1445;
elseif nargin==10
    MonitorOffset_x = Mon_Offset(1);
    MonitorOffset_z = Mon_Offset(2);
end


Monitor_xlim = [0 NumBlocks_col]*IlluminationBlock_size(1) + MonitorOffset_x;
Monitor_y = 0;
Monitor_zlim = [0 NumBlocks_row]*IlluminationBlock_size(2) + MonitorOffset_z;
N_monDiscr_row = Ndiscr_mon*NumBlocks_row;
N_monDiscr_col = Ndiscr_mon*NumBlocks_col;
Mon_xdiscr = linspace(Monitor_xlim(1),Monitor_xlim(2),N_monDiscr_col);
Mon_zdiscr = linspace(Monitor_zlim(2),Monitor_zlim(1),N_monDiscr_row);
BG_lightlevel = 0.00;


% SET MONITOR PATTERN
ActiveBlocks_matrix = zeros(NumBlocks_row,NumBlocks_col);
ActiveBlocks_matrix(ActiveBlocks_coord(:,1),ActiveBlocks_coord(:,2)) = 1; % Intensity of blocks.
MonitorPattern_comp = kron(ActiveBlocks_matrix,ones(Ndiscr_mon,Ndiscr_mon));
MonPattern = MonitorPattern_comp;
MonitorPattern_comp = (fliplr(MonitorPattern_comp))';
MonitorPattern_comp = MonitorPattern_comp(:)+BG_lightlevel;


[XX_wall,YY_wall,ZZ_wall] = meshgrid(Wall_xdiscr,LambertianWall_y,Wall_zdiscr);
[XX_mon,YY_mon,ZZ_mon] = meshgrid(Mon_xdiscr,Monitor_y,Mon_zdiscr);
Pos_wall = [XX_wall(:),YY_wall(:),ZZ_wall(:)];
Pos_mon = [XX_mon(:),YY_mon(:),ZZ_mon(:)];


Normal_monitor = [0, 1, 0];
Normal_wall = [0, -1, 0];

tempOnes = ones(prod(numPixels),1);
LL = zeros(prod(numPixels),1);

for ii=1:N_monDiscr_col*N_monDiscr_row
    if MonitorPattern_comp(ii)>0
        tDiff =(tempOnes*Pos_mon(ii,:) - Pos_wall);
        temp1 = (tDiff*Normal_wall')./(sum(tDiff.^2,2));
        temp1(temp1<0)= 0;
        temp2 = (-tDiff*Normal_monitor')./(sum(tDiff.^2,2).^2); %
        temp2(temp2<0)= 0;
        
        if Occ_present==true
            [Vx] = ComputeOccluderShadow(Pos_mon(ii,1),Wall_xdiscr,Dist_MonWall,Occluder(1,2),Occluder(:,1));
            [Vz] = ComputeOccluderShadow(Pos_mon(ii,3),Wall_zdiscr,Dist_MonWall,Occluder(2,2),Occluder(:,3));
            visibilityMat = (1-(1-Vz)'*(1-Vx))';
            visibilityMat = visibilityMat(:);
            if occ_stickpresent==true
                [Vx_s] = ComputeOccluderShadow(Pos_mon(ii,1),Wall_xdiscr,Dist_MonWall,occstick(1,2),occstick(:,1));
                [Vz_s] = ComputeOccluderShadow(Pos_mon(ii,3),Wall_zdiscr,Dist_MonWall,occstick(2,2),occstick(:,3));
                visibilityMat = (1- (((1-Vz)'*(1-Vx))|((1-Vz_s)'*(1-Vx_s))))';
                visibilityMat = visibilityMat(:);
            end
        else
            visibilityMat=1;
        end
        if viewAngleCorrection==true
            MM = ViewingAngleFactor(Pos_mon(ii,:),FOV,Dist_MonWall,numPixels);
        else
            MM=1;
        end
        LL(:,1) = LL(:,1) + temp1.*temp2.*visibilityMat*MonitorPattern_comp(ii).*(MM(:));
    end
end

Image = reshape(sum(LL,2),numPixels(1),numPixels(2)).';
Wall_discr.xdiscr = Wall_xdiscr;
Wall_discr.zdiscr = Wall_zdiscr;

end


function [M] = ViewingAngleFactor(MonitorPixel_xyz, Camera_FOV, D, numPixels)
%ViewingAngleFactor models the effect of the perceived brightness of the
% monitor changing depending on the angle from which it is viewed.

powexp = 2;20;%20; %5;
FOV_zdiscr = linspace(Camera_FOV(2,2),Camera_FOV(1,2),numPixels(2))';
FOV_xdiscr = linspace(Camera_FOV(2,1),Camera_FOV(1,1),numPixels(1))';
Mz = cos(1*atan((MonitorPixel_xyz(3)-FOV_zdiscr)./(D))).^powexp;
Mx = cos(atan((MonitorPixel_xyz(1)-FOV_xdiscr)./(D))).^0;%ones(size(FOV_xdiscr)); %cos(atan((MonitorPixel_xyz(1)-FOV_xdiscr)./(D))).^2;
% M = repmat(M,1,numPixels)';
M = (Mx*Mz');
end


function [Vx,X1,X2] = ComputeOccluderShadow(x_locs,xx,D_wall,Occl_d,Occ_edges)
%ComputeOccluderShadow 
%
%   Usage:
%       [Vx,X1,X2] = ComputeOccluderVisibility(x_locs,xx,D_wall,Occl_d,Occ_edges)
%
%   Input:
%       * x_locs:		.
%       * xx:           Dimension of the room
%       * D_wall:       Number of patches to divide wall (discr. of albedo).
%       * Occ_loc:      Number of illumation points (scalar).
%       * Occ_edges:    Number of detector points (scalar).
%   Output:
%       * Vx:           .

% Last Modified by $John Murray-Bruce$ at Boston University
% v1.0 28-Aug-2017 11:13:32


xx = xx(:).';              % Enforce a row vector.
N = length(xx);
[Ns] = length(x_locs);
XX = repmat(xx,Ns,1);

X1 = (D_wall*(Occ_edges(1) - x_locs)./Occl_d) +  x_locs;
X1 = repmat(X1,1,N);
X2 = (D_wall*(Occ_edges(2) - x_locs)./Occl_d) +  x_locs;
X2 = repmat(X2,1,N);

Vx = ones(Ns,N);
Vx( (XX>=X1)&(XX<=X2) ) =0;

end
